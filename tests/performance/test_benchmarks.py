"""Performance benchmarks for PyMatLib."""

import pytest
import time
import numpy as np
import sympy as sp
import os

from materforge.algorithms.interpolation import interpolate_value
from materforge.algorithms.piecewise_builder import PiecewiseBuilder

# Try to import psutil, skip memory tests if not available
try:
    import psutil
    PSUTIL_AVAILABLE = True
except ImportError:
    PSUTIL_AVAILABLE = False

def get_performance_threshold(base_threshold: float) -> float:
    """Adjust thresholds based on environment"""
    if os.getenv('CI'):
        # Allow 2x slower performance in CI
        return base_threshold * 2
    return base_threshold

class TestPerformance:
    """Performance benchmark tests."""
    @pytest.mark.slow
    def test_interpolation_performance(self):
        """Test interpolation performance with large datasets."""
        # Large dataset
        x_array = np.linspace(300, 1000, 10000)
        y_array = np.sin(x_array / 100) * 1000 + 1000
        # Time multiple interpolations
        start_time = time.time()
        for _ in range(1000):
            result = interpolate_value(500.0, x_array, y_array, 'constant', 'constant')
        end_time = time.time()
        avg_time = (end_time - start_time) / 1000
        # Adjust threshold for CI environment
        threshold = get_performance_threshold(0.005)  # 10ms in CI vs 5ms locally
        assert avg_time < threshold, f"Interpolation too slow: {avg_time:.6f}s per call (threshold: {threshold:.6f}s)"

    @pytest.mark.slow
    def test_piecewise_creation_performance(self):
        """Test piecewise function creation performance."""
        T = sp.Symbol('T')
        temp_array = np.linspace(300, 1000, 1000)
        prop_array = np.sin(temp_array / 100) * 1000 + 1000
        config = {'bounds': ['constant', 'constant']}
        start_time = time.time()
        piecewise_func = PiecewiseBuilder.build_from_data(
            temp_array, prop_array, T, config, "test_property"
        )
        end_time = time.time()
        creation_time = end_time - start_time
        # Adjust threshold for CI environment
        threshold = get_performance_threshold(2.0)  # 4s in CI vs 2s locally
        assert creation_time < threshold, f"Piecewise creation too slow: {creation_time:.3f}s (threshold: {threshold:.1f}s)"
        # Test evaluation performance
        start_time = time.time()
        for temp in [350, 450, 550, 650, 750]:
            result = float(piecewise_func.subs(T, temp))
        end_time = time.time()
        eval_time = (end_time - start_time) / 5
        # Adjust threshold for CI environment
        eval_threshold = get_performance_threshold(0.06)  # 120ms in CI vs 60ms locally
        assert eval_time < eval_threshold, f"Piecewise evaluation too slow: {eval_time:.6f}s per call (threshold: {eval_threshold:.6f}s)"

    @pytest.mark.slow
    def test_large_dataset_performance(self):
        """Test performance with very large datasets."""
        T = sp.Symbol('T')
        # Very large dataset
        temp_array = np.linspace(300, 2000, 5000)
        prop_array = np.random.normal(1000, 100, 5000)  # Random data
        config = {'bounds': ['constant', 'constant']}
        start_time = time.time()
        piecewise_func = PiecewiseBuilder.build_from_data(
            temp_array, prop_array, T, config, "large_property"
        )
        end_time = time.time()
        creation_time = end_time - start_time
        # Adjust threshold for CI environment
        threshold = get_performance_threshold(10.0)  # 20s in CI vs 10s locally
        assert creation_time < threshold, f"Large dataset processing too slow: {creation_time:.3f}s (threshold: {threshold:.1f}s)"
        # Test that the function is evaluable
        result = float(piecewise_func.subs(T, 1000))
        assert isinstance(result, float)
        assert not np.isnan(result)

    @pytest.mark.skipif(not PSUTIL_AVAILABLE, reason="psutil not available")
    def test_memory_efficiency(self):
        """Test memory usage doesn't grow excessively."""
        import os
        process = psutil.Process(os.getpid())
        initial_memory = process.memory_info().rss / 1024 / 1024  # MB
        T = sp.Symbol('T')
        # Create multiple piecewise functions
        for i in range(10):
            temp_array = np.linspace(300, 1000, 100)
            prop_array = np.random.normal(1000, 100, 100)
            config = {'bounds': ['constant', 'constant']}
            piecewise_func = PiecewiseBuilder.build_from_data(
                temp_array, prop_array, T, config, f"property_{i}"
            )
        final_memory = process.memory_info().rss / 1024 / 1024  # MB
        memory_increase = final_memory - initial_memory
        # Adjust threshold for CI environment
        threshold = get_performance_threshold(100)  # 200MB in CI vs 100MB locally
        assert memory_increase < threshold, f"Memory usage too high: {memory_increase:.2f}MB increase (threshold: {threshold:.0f}MB)"

    @pytest.mark.slow
    def test_material_creation_performance(self):
        """Test material creation performance."""
        import tempfile
        from pathlib import Path
        from ruamel.yaml import YAML
        from materforge.parsing.api import create_material
        T = sp.Symbol('T')
        # Complex material configuration
        config = {
            'name': 'Performance Test Steel',
            'material_type': 'alloy',
            'composition': {'Fe': 0.7, 'C': 0.1, 'Cr': 0.1, 'Ni': 0.1},
            'solidus_temperature': 1400.0,
            'liquidus_temperature': 1500.0,
            'initial_boiling_temperature': 2800.0,
            'final_boiling_temperature': 2900.0,
            'properties': {
                'density': 7850.0,
                'heat_capacity': {
                    'dependency': list(range(300, 1501, 50)),  # Many points
                    'value': [500 + i*0.1 for i in range(300, 1501, 50)],
                    'bounds': ['constant', 'constant']
                }
            }
        }
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
            yaml = YAML()
            yaml.dump(config, f)
            yaml_path = Path(f.name)
        try:
            start_time = time.time()
            material = create_material(yaml_path, T, enable_plotting=False)
            end_time = time.time()
            creation_time = end_time - start_time
            # Adjust threshold for CI environment
            threshold = get_performance_threshold(5.0)  # 10s in CI vs 5s locally
            assert creation_time < threshold, f"Material creation too slow: {creation_time:.3f}s (threshold: {threshold:.1f}s)"
            # Verify the material was created correctly
            assert material.name == "Performance Test Steel"
            assert hasattr(material, 'heat_capacity')
        finally:
            yaml_path.unlink()
