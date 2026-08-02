"""
Tests for config loading.

The interesting cases are the malformed ones. `configs/diseases/` ships five
empty files, so "empty YAML" is not a hypothetical input here.
"""

from pathlib import Path

import pytest

from syndrumnet.utils.config import Config, load_config, merge_configs

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_nested_access_both_ways(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("propagation:\n  alpha: 0.5\n  tolerance: 1.0e-6\nn_cores: 4\n")

    config = load_config(path)

    assert config.propagation.alpha == 0.5
    assert config['propagation']['tolerance'] == pytest.approx(1e-6)
    assert config.n_cores == 4


def test_empty_file_is_rejected(tmp_path):
    """
    yaml.safe_load returns None for an empty file. That used to reach Config
    and fail there with "'NoneType' object has no attribute 'items'", which
    says nothing about which file was at fault.
    """
    path = tmp_path / "empty.yaml"
    path.write_text("")

    with pytest.raises(ValueError, match="empty"):
        load_config(path)


def test_comments_only_file_is_rejected(tmp_path):
    path = tmp_path / "comments.yaml"
    path.write_text("# nothing but a comment\n")

    with pytest.raises(ValueError, match="empty"):
        load_config(path)


def test_non_mapping_file_is_rejected(tmp_path):
    path = tmp_path / "list.yaml"
    path.write_text("- one\n- two\n")

    with pytest.raises(ValueError, match="mapping"):
        load_config(path)


def test_missing_file_is_rejected(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_config(tmp_path / "does_not_exist.yaml")


def test_merge_configs_handles_dot_notation(tmp_path):
    path = tmp_path / "c.yaml"
    path.write_text("propagation:\n  alpha: 0.5\nn_cores: 4\n")

    merged = merge_configs(load_config(path), {'propagation.alpha': 0.7, 'n_cores': 8})

    assert merged.propagation.alpha == 0.7
    assert merged.n_cores == 8


def test_default_config_carries_the_visualization_block():
    """
    `scripts/make_figures.py` and `scripts/evaluate.py` read dpi, format and
    top_k from here. They were declared and read by nothing until now, so
    this guards against them drifting back out of use.
    """
    config = load_config(REPO_ROOT / "configs" / "default.yaml")

    assert isinstance(config.visualization, Config)
    assert config.visualization.get('dpi') == 300
    assert config.visualization.get('figure_format') == 'png'
    assert config.visualization.get('top_k_predictions') == 20
