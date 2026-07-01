# tests/test_python313.py
import sys
import pytest

@pytest.mark.skipif(sys.version_info < (3, 13), reason="Python 3.13+ only")
def test_python313_features():
    """Test Python 3.13 specific functionality."""
    import pastaq as pq
    # Your tests here
    assert hasattr(pq, '__version__')