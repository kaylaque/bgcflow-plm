from pathlib import Path

import pytest

FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir():
    return FIXTURES_DIR


@pytest.fixture
def mini_region_dir(tmp_path):
    """A tmp dir containing only mini_region.gbk (simulates a strain's antiSMASH dir)."""
    import shutil

    d = tmp_path / "antismash_strain"
    d.mkdir()
    shutil.copy(FIXTURES_DIR / "mini_region.gbk", d / "NODE_1.region001.gbk")
    return d


@pytest.fixture
def empty_region_dir(tmp_path):
    import shutil

    d = tmp_path / "antismash_empty"
    d.mkdir()
    shutil.copy(FIXTURES_DIR / "empty_region.gbk", d / "NODE_2.region001.gbk")
    return d


@pytest.fixture
def malformed_region_dir(tmp_path):
    import shutil

    d = tmp_path / "antismash_malformed"
    d.mkdir()
    shutil.copy(FIXTURES_DIR / "malformed_region.gbk", d / "NODE_3.region001.gbk")
    return d


@pytest.fixture
def mini_mibig_dir():
    return FIXTURES_DIR / "mini_mibig"
