"""Unit tests for property handlers."""

import sympy as sp
import tempfile
from pathlib import Path
from materforge.parsing.processors.property_handlers import (
    ConstantValuePropertyHandler,
    StepFunctionPropertyHandler,
    FileImportPropertyHandler,
    TabularDataPropertyHandler,
    PiecewiseEquationPropertyHandler
)

class TestConstantPropertyHandler:
    """Test cases for ConstantPropertyHandler."""
    def test_process_constant_property_float(self, sample_aluminum_element, temp_symbol):
        """Test processing constant float property."""
        from materforge.core.materials import Material

        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        handler = ConstantValuePropertyHandler()
        handler.set_processing_context(Path("."), None, set())
        handler.process_property(material, "test_density", 2700.0, temp_symbol)
        assert hasattr(material, "test_density")
        assert float(material.test_density) == 2700.0

    def test_process_constant_property_string(self, sample_aluminum_element, temp_symbol):
        """Test processing constant string property."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        handler = ConstantValuePropertyHandler()
        handler.set_processing_context(Path("."), None, set())
        handler.process_property(material, "test_property", "3.14", temp_symbol)
        assert hasattr(material, "test_property")
        assert float(material.test_property) == 3.14

class TestKeyValPropertyHandler:
    """Test cases for KeyValPropertyHandler."""
    def test_process_keyval_property(self, sample_aluminum_element, temp_symbol):
        """Test processing key-value property."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        config = {
            'dependency': [300, 400, 500],
            'value': [900, 950, 1000],
            'bounds': ['constant', 'constant']
        }
        handler = TabularDataPropertyHandler()
        handler.set_processing_context(Path("."), None, set())
        handler.process_property(material, "heat_capacity", config, temp_symbol)
        assert hasattr(material, "heat_capacity")
        assert isinstance(material.heat_capacity, sp.Piecewise)
        # Test evaluation
        result = float(material.heat_capacity.subs(temp_symbol, 350))
        assert 900 < result < 1000  # Should be interpolated

class TestStepFunctionPropertyHandler:
    """Test cases for StepFunctionPropertyHandler."""
    def test_process_step_function_property(self, sample_aluminum_element, temp_symbol):
        """Test processing step function property."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        config = {
            'dependency': 'melting_temperature',
            'value': [100.0, 200.0]
        }
        handler = StepFunctionPropertyHandler()
        handler.set_processing_context(Path("."), None, set())
        handler.process_property(material, "test_step", config, temp_symbol)
        assert hasattr(material, "test_step")
        assert isinstance(material.test_step, sp.Piecewise)
        # Test evaluation below and above transition
        result_below = float(material.test_step.subs(temp_symbol, 800))  # Below melting point
        result_above = float(material.test_step.subs(temp_symbol, 1200))  # Above melting point
        assert result_below == 100.0
        assert result_above == 200.0

class TestPiecewiseEquationPropertyHandler:
    """Test cases for PiecewiseEquationPropertyHandler."""
    def test_process_piecewise_equation_property(self, sample_aluminum_element, temp_symbol):
        """Test processing piecewise equation property."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        config = {
            'dependency': [300, 500, 700],
            'equation': ['2*T + 100', '3*T - 200'],
            'bounds': ['constant', 'constant']
        }
        handler = PiecewiseEquationPropertyHandler()
        handler.set_processing_context(Path("."), None, set())
        handler.process_property(material, "test_piecewise", config, temp_symbol)
        assert hasattr(material, "test_piecewise")
        assert isinstance(material.test_piecewise, sp.Piecewise)
        # Test evaluation in different segments
        result_first = float(material.test_piecewise.subs(temp_symbol, 400))  # First segment
        result_second = float(material.test_piecewise.subs(temp_symbol, 600))  # Second segment
        expected_first = 2*400 + 100  # 900
        expected_second = 3*600 - 200  # 1600
        assert abs(result_first - expected_first) < 1e-10
        assert abs(result_second - expected_second) < 1e-10

class TestFilePropertyHandler:
    """Test cases for FilePropertyHandler."""
    def test_process_file_property_csv(self, sample_aluminum_element, temp_symbol):
        """Test processing property from CSV file."""
        from materforge.core.materials import Material
        material = Material(
            name="Test Material",
            material_type="pure_metal",
            elements=[sample_aluminum_element],
            composition=[1.0],
            melting_temperature=sp.Float(933.47),
            boiling_temperature=sp.Float(2792.0)
        )
        # Create temporary CSV file
        csv_content = """temperature,heat_capacity
                        300,900
                        400,950
                        500,1000
                        """
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            csv_path = Path(f.name)
        try:
            config = {
                'file_path': csv_path.name,
                'dependency_column': 'temperature',
                'property_column': 'heat_capacity',
                'bounds': ['constant', 'constant']
            }
            handler = FileImportPropertyHandler()
            handler.set_processing_context(csv_path.parent, None, set())
            handler.process_property(material, "file_heat_capacity", config, temp_symbol)
            assert hasattr(material, "file_heat_capacity")
            assert isinstance(material.file_heat_capacity, sp.Piecewise)
            # Test evaluation
            result = float(material.file_heat_capacity.subs(temp_symbol, 350))
            assert 900 < result < 1000  # Should be interpolated
        finally:
            csv_path.unlink()
