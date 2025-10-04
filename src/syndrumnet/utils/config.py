"""
Configuration management for SyndrumNET pipeline.

Handles loading YAML configs, merging defaults with overrides, and 
providing nested dictionary access.
"""

from pathlib import Path
from typing import Any, Dict, Optional, Union
import yaml


class Config:
    """
    Configuration container with nested access support.
    
    Supports both dict-like and attribute-like access:
        config['data']['network_sources']  # dict-style
        config.data.network_sources         # attribute-style
    """
    
    def __init__(self, config_dict: Dict[str, Any]) -> None:
        """Initialize config from dictionary."""
        self._config = config_dict
        
        # Convert nested dicts to Config objects
        for key, value in config_dict.items():
            if isinstance(value, dict):
                setattr(self, key, Config(value))
            else:
                setattr(self, key, value)
    
    def __getitem__(self, key: str) -> Any:
        """Dict-like access."""
        return self._config[key]
    
    def __setitem__(self, key: str, value: Any) -> None:
        """Dict-like assignment."""
        self._config[key] = value
        setattr(self, key, value)
    
    def get(self, key: str, default: Any = None) -> Any:
        """Get with default value."""
        return self._config.get(key, default)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert back to plain dictionary."""
        result = {}
        for key, value in self._config.items():
            if isinstance(value, Config):
                result[key] = value.to_dict()
            else:
                result[key] = value
        return result
    
    def __repr__(self) -> str:
        return f"Config({self._config})"


def load_config(config_path: Union[str, Path]) -> Config:
    """
    Load configuration from YAML file.
    
    Parameters
    ----------
    config_path : str or Path
        Path to YAML configuration file.
        
    Returns
    -------
    Config
        Configuration object with nested access.
        
    Examples
    --------
    >>> config = load_config('configs/default.yaml')
    >>> config.propagation.alpha
    0.5
    >>> config['propagation']['max_iterations']
    1000
    """
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config_dict = yaml.safe_load(f)
    
    return Config(config_dict)


def merge_configs(base: Config, override: Dict[str, Any]) -> Config:
    """
    Merge override dictionary into base config.
    
    Supports nested keys with dot notation: 'propagation.alpha' = 0.7
    
    Parameters
    ----------
    base : Config
        Base configuration.
    override : dict
        Dictionary of override values. Keys can use dot notation.
        
    Returns
    -------
    Config
        Merged configuration.
        
    Examples
    --------
    >>> config = load_config('default.yaml')
    >>> overrides = {'propagation.alpha': 0.7, 'n_cores': 8}
    >>> config = merge_configs(config, overrides)
    """
    merged = base.to_dict()
    
    for key, value in override.items():
        # Handle dot notation: 'propagation.alpha' -> ['propagation', 'alpha']
        keys = key.split('.')
        target = merged
        
        for k in keys[:-1]:
            if k not in target:
                target[k] = {}
            target = target[k]
        
        target[keys[-1]] = value
    
    return Config(merged)


def save_config(config: Config, output_path: Union[str, Path]) -> None:
    """
    Save configuration to YAML file.
    
    Parameters
    ----------
    config : Config
        Configuration to save.
    output_path : str or Path
        Output YAML file path.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        yaml.dump(config.to_dict(), f, default_flow_style=False, sort_keys=False)