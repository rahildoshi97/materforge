"""Tests for data handler module with improved error handling."""

import pytest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from materforge.parsing.io.data_handler import (
    load_property_data, _extract_data_columns, _clean_and_validate_data
)


class TestLoadPropertyData:
    """Test the main data loading function."""

    def test_load_csv_data(self):
        """Test loading data from CSV file."""
        csv_content = """Temperature,Property
300,100
400,150
500,200"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            csv_path = Path(f.name)
        try:
            config = {
                'file_path': str(csv_path),
                'dependency_column': 'Temperature',
                'property_column': 'Property'
            }
            temp_array, prop_array = load_property_data(config)
            assert len(temp_array) == 3
            assert len(prop_array) == 3
            assert temp_array[0] == 300
            assert prop_array[0] == 100
        finally:
            csv_path.unlink()

    def test_load_excel_data(self):
        """Test loading data from Excel file."""
        df = pd.DataFrame({
            'Temp': [300, 400, 500],
            'Value': [100, 150, 200]
        })
        with tempfile.NamedTemporaryFile(suffix='.xlsx', delete=False) as f:
            excel_path = Path(f.name)
        try:
            df.to_excel(excel_path, index=False)
            config = {
                'file_path': str(excel_path),
                'dependency_column': 'Temp',
                'property_column': 'Value'
            }
            temp_array, prop_array = load_property_data(config)
            assert len(temp_array) == 3
            assert len(prop_array) == 3
        finally:
            excel_path.unlink()

    def test_load_data_missing_file(self):
        """Test error handling for missing files."""
        config = {
            'file_path': 'nonexistent.csv',
            'dependency_column': 'Temperature',
            'property_column': 'Property'
        }
        with pytest.raises(FileNotFoundError):
            load_property_data(config)

    def test_load_data_unsupported_format(self):
        """Test error handling for unsupported file formats."""
        with tempfile.NamedTemporaryFile(suffix='.xyz', delete=False) as f:  # Use unsupported extension
            xyz_path = Path(f.name)
        try:
            config = {
                'file_path': str(xyz_path),
                'dependency_column': 'Temperature',
                'property_column': 'Property'
            }
            with pytest.raises(ValueError, match="Unsupported file type"):
                load_property_data(config)
        finally:
            xyz_path.unlink()

    def test_load_data_with_missing_values(self):
        """Test handling of missing values in data."""
        csv_content = """Temperature,Property
300,100
400,
500,200
,150"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            csv_path = Path(f.name)
        try:
            config = {
                'file_path': str(csv_path),
                'dependency_column': 'Temperature',
                'property_column': 'Property'
            }
            temp_array, prop_array = load_property_data(config)
            # Should have cleaned out rows with missing values
            assert len(temp_array) == 2
            assert len(prop_array) == 2
            assert temp_array[0] == 300
            assert prop_array[0] == 100
        finally:
            csv_path.unlink()

    def test_load_data_column_not_found(self):
        """Test error handling for missing columns."""
        csv_content = """Temperature,Property
300,100
400,150"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            csv_path = Path(f.name)
        try:
            config = {
                'file_path': str(csv_path),
                'dependency_column': 'NonexistentColumn',
                'property_column': 'Property'
            }
            with pytest.raises(ValueError, match="not found"):
                load_property_data(config)
        finally:
            csv_path.unlink()

    def test_load_empty_text_file(self):
        """Test handling of empty text files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            # Create empty file
            pass
        try:
            config = {
                'file_path': f.name,
                'dependency_column': 'Temperature',
                'property_column': 'Property'
            }
            with pytest.raises(ValueError, match="empty|No data"):
                load_property_data(config)
        finally:
            Path(f.name).unlink()


class TestDataValidation:
    """Test data validation functions."""

    def test_clean_and_validate_data_valid(self):
        """Test cleaning valid data."""
        temp = np.array([300, 400, 500])
        prop = np.array([100, 150, 200])
        clean_temp, clean_prop = _clean_and_validate_data(temp, prop, "test.csv")
        assert len(clean_temp) == 3
        assert len(clean_prop) == 3
        np.testing.assert_array_equal(clean_temp, temp)
        np.testing.assert_array_equal(clean_prop, prop)

    def test_clean_and_validate_data_with_acceptable_nans(self):
        """Test cleaning data with acceptable NaN values (under 50% threshold)."""
        temp = np.array([300, np.nan, 500, 600])  # 25% missing
        prop = np.array([100, 150, np.nan, 200])  # 25% missing
        clean_temp, clean_prop = _clean_and_validate_data(temp, prop, "test.csv")
        assert len(clean_temp) == 2  # Only 2 valid rows remain
        assert len(clean_prop) == 2
        assert clean_temp[0] == 300
        assert clean_prop[0] == 100

    def test_clean_and_validate_data_empty(self):
        """Test error handling for empty data."""
        temp = np.array([])
        prop = np.array([])
        with pytest.raises(ValueError, match="No valid data found"):
            _clean_and_validate_data(temp, prop, "test.csv")

    def test_clean_and_validate_data_excessive_missing(self):
        """Test error handling for excessive missing values."""
        temp = np.array([300, np.nan, np.nan])  # 66.7% missing
        prop = np.array([100, np.nan, np.nan])  # 66.7% missing
        with pytest.raises(ValueError, match="Too many missing values"):
            _clean_and_validate_data(temp, prop, "test.csv")

    def test_clean_and_validate_data_duplicates(self):
        """Test handling of duplicate temperature values."""
        temp = np.array([300, 300, 400, 500])
        prop = np.array([100, 105, 150, 200])
        clean_temp, clean_prop = _clean_and_validate_data(temp, prop, "test.csv")
        # Should remove duplicates, keeping first occurrence
        assert len(clean_temp) == 3
        assert clean_temp[0] == 300
        assert clean_prop[0] == 100  # First occurrence kept

    def test_clean_and_validate_insufficient_data_after_cleaning(self):
        """Test error when insufficient data remains after cleaning."""
        temp = np.array([300])  # Only one data point
        prop = np.array([100])
        with pytest.raises(ValueError, match="Insufficient valid data points"):
            _clean_and_validate_data(temp, prop, "test.csv")


class TestConfigValidation:
    """Test configuration validation."""

    def test_missing_config_keys(self):
        """Test error handling for missing configuration keys."""
        config = {
            'file_path': 'test.csv',
            # Missing dependency_column and property_column
        }
        with pytest.raises(ValueError, match="Missing required configuration keys"):
            load_property_data(config)

    def test_empty_file_path(self):
        """Test error handling for empty file path."""
        config = {
            'file_path': '',
            'dependency_column': 'Temperature',
            'property_column': 'Property'
        }
        with pytest.raises(ValueError, match="File path cannot be empty"):
            load_property_data(config)
