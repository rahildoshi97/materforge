import logging
import numpy as np
from typing import Union, Tuple, Dict
import pandas as pd
from pathlib import Path

from materforge.parsing.config.yaml_keys import FILE_PATH_KEY, DEPENDENCY_COLUMN_KEY, PROPERTY_COLUMN_KEY
from materforge.data.constants import ProcessingConstants, FileConstants

logger = logging.getLogger(__name__)


def load_property_data(file_config: Dict[str, Union[str, int]], header: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Reads temperature and property data from a file with comprehensive error handling.
    Args:
        file_config: Dictionary containing file configuration with keys:
            - file_path: Path to data file
            - temperature_header: Temperature column name/index
            - value_header: Property column name/index
        header: Indicates if the file contains a header row
    Returns:
        Tuple of (temperature_array, property_array) as numpy arrays
    Raises:
        FileNotFoundError: If the specified file doesn't exist
        ValueError: If data validation fails or file format is unsupported
        PermissionError: If file cannot be read due to permissions
    """
    # Validate input configuration
    _validate_file_config(file_config)
    # Extract configuration
    file_path = Path(file_config[FILE_PATH_KEY])
    temp_col = file_config[DEPENDENCY_COLUMN_KEY]
    prop_col = file_config[PROPERTY_COLUMN_KEY]
    # Check file existence and permissions
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    if not file_path.is_file():
        raise ValueError(f"Path is not a file: {file_path}")
    # Validate file extension early
    file_extension = file_path.suffix.lower()
    if file_extension not in FileConstants.SUPPORTED_EXTENSIONS:
        raise ValueError(f"Unsupported file type: '{file_extension}'. "
                         f"Supported types are: {FileConstants.SUPPORTED_EXTENSIONS}")
    # Check file size
    file_size_mb = file_path.stat().st_size / (1024 * 1024)
    if file_size_mb > FileConstants.MAX_FILE_SIZE_MB:
        raise ValueError(f"File size ({file_size_mb:.2f} MB) exceeds the maximum limit "
                         f"of {FileConstants.MAX_FILE_SIZE_MB} MB.")
    try:
        # Read file based on extension
        if file_extension == '.xlsx':
            df = _read_excel_file(file_path, header, temp_col, prop_col)
        elif file_extension == '.csv':
            df = _read_csv_file(file_path, header, temp_col, prop_col)
        else:  # .txt files
            return _read_text_file(file_path, header, temp_col, prop_col)
    except PermissionError as e:
        raise PermissionError(f"Permission denied reading file {file_path}: {str(e)}") from e
    except Exception as e:
        raise ValueError(f"Error reading file {file_path}: {str(e)}") from e
    # Extract and validate data
    temp_array, prop_array = _extract_data_columns(df, temp_col, prop_col, str(file_path))
    temp_array, prop_array = _clean_and_validate_data(temp_array, prop_array, str(file_path))
    return temp_array, prop_array


def _validate_file_config(file_config: Dict) -> None:
    """Validate the file configuration dictionary."""
    required_keys = {FILE_PATH_KEY, DEPENDENCY_COLUMN_KEY, PROPERTY_COLUMN_KEY}
    missing_keys = required_keys - set(file_config.keys())
    if missing_keys:
        raise ValueError(f"Missing required configuration keys: {missing_keys}")
    if not file_config[FILE_PATH_KEY]:
        raise ValueError("File path cannot be empty")


def _read_excel_file(file_path: Path, header: bool, temp_col: Union[str, int],
                     prop_col: Union[str, int]) -> pd.DataFrame:
    """Read Excel file with proper error handling."""
    try:
        return pd.read_excel(
            file_path,
            header=0 if header else None,
            na_values=FileConstants.NA_VALUES,
            converters={
                temp_col: lambda x: pd.to_numeric(x, errors='coerce'),
                prop_col: lambda x: pd.to_numeric(x, errors='coerce')
            } if header else None
        )
    except ImportError as e:
        raise ValueError("Excel file support requires openpyxl. Install with: pip install openpyxl") from e


def _read_csv_file(file_path: Path, header: bool, temp_col: Union[str, int],
                   prop_col: Union[str, int]) -> pd.DataFrame:
    """Read CSV file with proper error handling."""
    return pd.read_csv(
        file_path,
        header=0 if header else None,
        na_values=FileConstants.NA_VALUES,
        encoding=FileConstants.DEFAULT_ENCODING,
        converters={
            temp_col: lambda x: pd.to_numeric(x, errors='coerce'),
            prop_col: lambda x: pd.to_numeric(x, errors='coerce')
        } if header else None
    )


def _read_text_file(file_path: Path, header: bool, temp_col: Union[str, int],
                    prop_col: Union[str, int]) -> Tuple[np.ndarray, np.ndarray]:
    """Read text file with comprehensive error handling."""
    try:
        if header:
            # Read header to get column names
            with open(file_path, 'r', encoding=FileConstants.DEFAULT_ENCODING) as f:
                header_line = f.readline().strip()
                if not header_line:
                    raise ValueError("File appears to be empty or has no header")
                column_names = header_line.split()
            # Read data
            data = pd.read_csv(
                file_path,
                sep=r'\s+',
                skiprows=1,
                header=None,
                na_values=FileConstants.NA_VALUES,
                encoding=FileConstants.DEFAULT_ENCODING,
                engine='python'  # Explicitly specify engine for regex separator
            ).values
            if data.size == 0:
                raise ValueError("No data found in file after header")
            # Handle column mapping
            temp_idx = _get_column_index(temp_col, column_names, data.shape[1], "temperature")
            prop_idx = _get_column_index(prop_col, column_names, data.shape[1], "property")
        else:
            # No header case
            data = pd.read_csv(
                file_path,
                sep=r'\s+',
                header=None,
                na_values=FileConstants.NA_VALUES,
                encoding=FileConstants.DEFAULT_ENCODING,
                engine='python'
            ).values
            if isinstance(temp_col, str) or isinstance(prop_col, str):
                raise ValueError("Column names specified, but file has no header row")
            temp_idx, prop_idx = temp_col, prop_col
        temp = data[:, temp_idx]
        prop = data[:, prop_idx]
        return temp, prop
    except pd.errors.EmptyDataError as e:
        raise ValueError(f"No data found in text file: {str(e)}") from e
    except Exception as e:
        raise ValueError(f"Error processing text file: {str(e)}") from e


def _extract_data_columns(df: pd.DataFrame, temp_col: Union[str, int], prop_col: Union[str, int],
                          file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Extract temperature and property columns from DataFrame with enhanced error handling."""
    if df.empty:
        raise ValueError(f"No data found in file: {file_path}")
    # Extract temperature column
    temp_series = _extract_column(df, temp_col, "temperature", file_path)
    # Extract property column
    prop_series = _extract_column(df, prop_col, "property", file_path)
    # Convert to numeric arrays
    temp_array, prop_array = _convert_to_numeric_arrays(temp_series, prop_series, file_path)
    return temp_array, prop_array


def _extract_column(df: pd.DataFrame, col_identifier: Union[str, int],
                    col_type: str, file_path: str) -> pd.Series:
    """Extract a single column from DataFrame with proper error handling."""
    if isinstance(col_identifier, str):
        if col_identifier not in df.columns:
            available_cols = ', '.join(df.columns.astype(str))
            raise ValueError(f"{col_type.capitalize()} column '{col_identifier}' not found in file {file_path}. "
                             f"Available columns: {available_cols}")
        return df[col_identifier]
    else:
        if col_identifier >= len(df.columns):
            raise ValueError(f"{col_type.capitalize()} column index {col_identifier} out of bounds "
                             f"(file has {len(df.columns)} columns)")
        return df.iloc[:, col_identifier]


def _convert_to_numeric_arrays(temp_series: pd.Series, prop_series: pd.Series,
                               file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Convert pandas Series to numeric numpy arrays with validation."""
    try:
        # Convert to numeric, handling any remaining non-numeric values
        temp_numeric = pd.to_numeric(temp_series, errors='coerce')
        prop_numeric = pd.to_numeric(prop_series, errors='coerce')
        # Convert to float64 numpy arrays
        temp_array = np.asarray(temp_numeric, dtype=np.float64)
        prop_array = np.asarray(prop_numeric, dtype=np.float64)
        # Validate conversion success
        if not np.issubdtype(temp_array.dtype, np.number):
            raise ValueError(f"Temperature column could not be converted to numeric type."
                             f"Got dtype: {temp_array.dtype}")
        if not np.issubdtype(prop_array.dtype, np.number):
            raise ValueError(f"Property column could not be converted to numeric type. Got dtype: {prop_array.dtype}")
        # Log conversion statistics
        temp_nan_count = np.sum(np.isnan(temp_array))
        prop_nan_count = np.sum(np.isnan(prop_array))
        if temp_nan_count > 0:
            logger.warning(f"Temperature column has {temp_nan_count} NaN values after conversion")
        if prop_nan_count > 0:
            logger.warning(f"Property column has {prop_nan_count} NaN values after conversion")
        return temp_array, prop_array
    except Exception as e:
        raise ValueError(f"Failed to convert data to numeric format in {file_path}: {str(e)}") from e


def _get_column_index(col_identifier: Union[str, int], column_names: list,
                      num_cols: int, col_type: str) -> int:
    """Get column index from name or validate numeric index."""
    if isinstance(col_identifier, str):
        try:
            return column_names.index(col_identifier)
        except ValueError:
            raise ValueError(f"{col_type.capitalize()} column '{col_identifier}' not found. "
                             f"Available columns: {', '.join(column_names)}")
    else:
        if col_identifier >= num_cols:
            raise ValueError(f"{col_type.capitalize()} column index {col_identifier} "
                             f"out of bounds (file has {num_cols} columns)")
        return col_identifier


def _clean_and_validate_data(temp_array: np.ndarray, prop_array: np.ndarray,
                             file_path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Validate data quality and handle missing values with improved logic."""
    # Check for completely empty arrays
    if len(temp_array) == 0 or len(prop_array) == 0:
        raise ValueError(f"No valid data found in file: {file_path}")
    # Identify rows with missing values
    temp_nan_mask = np.isnan(temp_array)
    prop_nan_mask = np.isnan(prop_array)
    any_nan_mask = temp_nan_mask | prop_nan_mask
    # Handle missing values if found
    if np.any(any_nan_mask):
        nan_count = np.sum(any_nan_mask)
        total_count = len(temp_array)
        nan_percentage = (nan_count / total_count) * 100
        logger.warning(f"Found {nan_count} rows ({nan_percentage:.1f}%) with missing values in {file_path}")
        # Check if too many missing values
        if nan_percentage > ProcessingConstants.MAX_MISSING_VALUE_PERCENTAGE:
            raise ValueError(f"Too many missing values ({nan_percentage:.1f}%) in file: {file_path}. "
                             "Please clean the data or check file format.")
        # Remove rows with missing values
        valid_mask = ~any_nan_mask
        temp_array = temp_array[valid_mask]
        prop_array = prop_array[valid_mask]
        logger.info(f"Removed {nan_count} rows with missing values. "
                    f"Remaining data points: {len(temp_array)}")
    # Final validation
    if len(temp_array) < ProcessingConstants.MIN_DATA_POINTS:
        raise ValueError(f"Insufficient valid data points ({len(temp_array)}) after cleaning missing values. "
                         f"Minimum required: {ProcessingConstants.MIN_DATA_POINTS}")
    # Handle duplicate temperatures
    temp_array, prop_array = _remove_duplicate_temperatures(temp_array, prop_array)
    # Ensure proper ordering
    temp_array, prop_array = _ensure_temperature_ordering(temp_array, prop_array)
    return temp_array, prop_array


def _remove_duplicate_temperatures(temp_array: np.ndarray, prop_array: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Remove duplicate temperature entries, keeping the first occurrence."""
    unique_temp, unique_indices = np.unique(temp_array, return_index=True)
    if len(unique_temp) < len(temp_array):
        duplicate_count = len(temp_array) - len(unique_temp)
        logger.warning(f"Found {duplicate_count} duplicate temperature entries. Removing duplicates.")
        # Sort indices to maintain original order
        unique_indices = np.sort(unique_indices)
        temp_array = temp_array[unique_indices]
        prop_array = prop_array[unique_indices]
    return temp_array, prop_array


def _ensure_temperature_ordering(temp_array: np.ndarray, prop_array: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Ensure temperature array is in ascending order."""
    if not np.all(np.diff(temp_array) > 0):
        logger.info("Sorting data by temperature")
        sort_indices = np.argsort(temp_array)
        temp_array = temp_array[sort_indices]
        prop_array = prop_array[sort_indices]
    return temp_array, prop_array
