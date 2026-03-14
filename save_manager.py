"""
Save Manager for NV Experiment Data and Parameters.

This module provides unified save/load functionality for:
- NumPy arrays (.npy) for data storage
- YAML (.yaml) for human-readable parameter files (separate from data)
- JSON (.json) as alternative parameter format
- Automatic folder creation with date/sequence structure

Date: December 2024
"""

import numpy as np
import yaml
import json
import time
from pathlib import Path
from typing import Union, Dict, Any, Optional


class SaveManager:
    """
    Manages saving and loading of experiment data and parameters.

    Key Features:
    - NumPy arrays for data (backward compatible)
    - Separate YAML/JSON files for parameters (human-readable)
    - Automatic folder creation: Saved_Data/YYYY-MM-DD/sequence_XXX/
    - NumPy type conversion for serialization

    File Structure:
        Saved_Data/
        └── 2024-12-15/
            ├── esr_001/
            │   ├── data_001.npy      # Raw data array
            │   └── params_001.yaml   # Parameters (separate!)
            └── rabi_002/
                ├── data_002.npy
                └── params_002.yaml
    """

    def __init__(self, base_dir: Optional[Path] = None):
        """
        Initialize SaveManager.

        Args:
            base_dir: Base directory for saved data. If None, uses ../Saved_Data
        """
        if base_dir is None:
            # Default: one level up from current directory
            self.base_dir = Path.cwd().parent / "Saved_Data"
        else:
            self.base_dir = Path(base_dir)

        # Ensure base directory exists
        self.base_dir.mkdir(parents=True, exist_ok=True)

    # ========================================================================
    # SAVE METHODS
    # ========================================================================

    def save_numpy(self, data: np.ndarray, folder_number: str, sequence: str):
        """
        Save numpy data array.

        Args:
            data: NumPy array to save (typically shape: (Nruns, Nscanpts, samples))
            folder_number: Folder number identifier (e.g., "001", "002")
            sequence: Sequence type (e.g., "esr", "rabi", "t2")

        File saved to: Saved_Data/YYYY-MM-DD/sequence_XXX/data_XXX.npy
        """
        file_dir = self._create_folder(sequence, folder_number)
        data_file = file_dir / f"data_{folder_number}.npy"

        np.save(data_file, data, allow_pickle=False)
        print(f"✔ Data saved to {data_file}")

    def save_params_yaml(self, params: dict, folder_number: str, sequence: str):
        """
        Save parameters as YAML (separate from data file).

        This creates a human-readable parameter file that can be viewed
        without loading the data array.

        Args:
            params: Parameter dictionary (typically params_dict from config)
            folder_number: Folder number identifier (e.g., "001", "002")
            sequence: Sequence type (e.g., "esr", "rabi", "t2")

        File saved to: Saved_Data/YYYY-MM-DD/sequence_XXX/params_XXX.yaml
        """
        file_dir = self._create_folder(sequence, folder_number)
        yaml_file = file_dir / f"params_{folder_number}.yaml"

        # Convert numpy types to Python types for serialization
        params_serializable = self._convert_numpy_to_python(params)

        with open(yaml_file, 'w') as f:
            yaml.dump(params_serializable, f, indent=4, default_flow_style=False)

        print(f"✔ Parameters saved to {yaml_file}")

    def save_params_json(self, params: dict, folder_number: str, sequence: str):
        """
        Save parameters as JSON (alternative to YAML).

        Args:
            params: Parameter dictionary
            folder_number: Folder number identifier (e.g., "001", "002")
            sequence: Sequence type (e.g., "esr", "rabi", "t2")

        File saved to: Saved_Data/YYYY-MM-DD/sequence_XXX/params_XXX.json
        """
        file_dir = self._create_folder(sequence, folder_number)
        json_file = file_dir / f"params_{folder_number}.json"

        # Convert numpy types to Python types for serialization
        params_serializable = self._convert_numpy_to_python(params)

        with open(json_file, 'w') as f:
            json.dump(params_serializable, f, indent=4)

        print(f"✔ Parameters saved to {json_file}")

    def save_experiment(self, data: np.ndarray, params: dict,
                       folder_number: str, sequence: str,
                       param_format: str = 'yaml'):
        """
        Save both data and parameters together.

        Convenience method that saves data as .npy and parameters as YAML/JSON.

        Args:
            data: NumPy array to save
            params: Parameter dictionary
            folder_number: Folder number identifier
            sequence: Sequence type
            param_format: 'yaml' or 'json' for parameter file format

        Example:
            save_manager.save_experiment(
                data=data_array,
                params=params_dict,
                folder_number="001",
                sequence="esr",
                param_format='yaml'
            )
        """
        self.save_numpy(data, folder_number, sequence)

        if param_format.lower() == 'yaml':
            self.save_params_yaml(params, folder_number, sequence)
        elif param_format.lower() == 'json':
            self.save_params_json(params, folder_number, sequence)
        else:
            raise ValueError(f"Invalid param_format: {param_format}. Use 'yaml' or 'json'.")

    # ========================================================================
    # LOAD METHODS
    # ========================================================================

    def load_numpy(self, filepath: Union[str, Path]) -> np.ndarray:
        """
        Load numpy data array.

        Args:
            filepath: Path to .npy file

        Returns:
            NumPy array
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Data file not found: {filepath}")

        data = np.load(filepath, allow_pickle=False)
        print(f"✔ Data loaded from {filepath}")
        return data

    def load_params(self, filepath: Union[str, Path]) -> dict:
        """
        Load parameters from YAML or JSON file.

        Automatically detects file format from extension.

        Args:
            filepath: Path to .yaml or .json file

        Returns:
            Parameter dictionary
        """
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Parameter file not found: {filepath}")

        if filepath.suffix == '.yaml' or filepath.suffix == '.yml':
            with open(filepath, 'r') as f:
                params = yaml.safe_load(f)
            print(f"✔ Parameters loaded from {filepath}")
            return params

        elif filepath.suffix == '.json':
            with open(filepath, 'r') as f:
                params = json.load(f)
            print(f"✔ Parameters loaded from {filepath}")
            return params

        else:
            raise ValueError(f"Unsupported file format: {filepath.suffix}. Use .yaml or .json")

    def load_experiment(self, data_filepath: Union[str, Path],
                       params_filepath: Optional[Union[str, Path]] = None) -> tuple:
        """
        Load both data and parameters together.

        If params_filepath is not provided, attempts to find it automatically
        by replacing data file extension with .yaml or .json.

        Args:
            data_filepath: Path to .npy data file
            params_filepath: Optional path to parameter file. If None, auto-detected.

        Returns:
            Tuple of (data_array, params_dict)

        Example:
            data, params = save_manager.load_experiment(
                "Saved_Data/2024-12-15/esr_001/data_001.npy"
            )
            # Auto-loads params_001.yaml from same folder
        """
        data_filepath = Path(data_filepath)
        data = self.load_numpy(data_filepath)

        # Auto-detect parameter file if not provided
        if params_filepath is None:
            # Try YAML first, then JSON
            yaml_file = data_filepath.parent / data_filepath.name.replace('.npy', '.yaml')
            json_file = data_filepath.parent / data_filepath.name.replace('.npy', '.json')

            if yaml_file.exists():
                params_filepath = yaml_file
            elif json_file.exists():
                params_filepath = json_file
            else:
                raise FileNotFoundError(
                    f"Could not find parameter file (tried {yaml_file.name} and {json_file.name})"
                )

        params = self.load_params(params_filepath)

        return data, params

    # ========================================================================
    # HELPER METHODS
    # ========================================================================

    def _create_folder(self, sequence: str, folder_number: str) -> Path:
        """
        Create save directory with date/sequence structure.

        Structure: Saved_Data/YYYY-MM-DD/sequence_XXX/

        Args:
            sequence: Sequence type (e.g., "esr", "rabi")
            folder_number: Folder number (e.g., "001")

        Returns:
            Path to created folder
        """
        # Create date directory (e.g., 2024-12-15)
        date_dir = self.base_dir / time.strftime("%Y-%m-%d", time.localtime())
        date_dir.mkdir(parents=True, exist_ok=True)

        # Create sequence folder (e.g., esr_001)
        file_dir = date_dir / f"{sequence}_{folder_number}"
        file_dir.mkdir(exist_ok=True)

        return file_dir

    def _convert_numpy_to_python(self, obj: Any) -> Any:
        """
        Recursively convert numpy types to Python types for serialization.

        Handles:
        - np.ndarray → list
        - np.int64, np.float64, etc. → int, float
        - dict → recursively converted dict
        - list → recursively converted list
        - Path → str

        Args:
            obj: Object to convert

        Returns:
            Converted object with Python-native types
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.generic):
            return obj.item()
        elif isinstance(obj, dict):
            return {key: self._convert_numpy_to_python(value)
                    for key, value in obj.items()}
        elif isinstance(obj, list):
            return [self._convert_numpy_to_python(item) for item in obj]
        elif isinstance(obj, Path):
            return str(obj)
        else:
            return obj

    # ========================================================================
    # UTILITY METHODS
    # ========================================================================

    def list_experiments(self, date: Optional[str] = None) -> list:
        """
        List all experiment folders.

        Args:
            date: Optional date string (YYYY-MM-DD). If None, uses today's date.

        Returns:
            List of experiment folder paths

        Example:
            experiments = save_manager.list_experiments("2024-12-15")
            # Returns: [esr_001, esr_002, rabi_001, ...]
        """
        if date is None:
            date = time.strftime("%Y-%m-%d", time.localtime())

        date_dir = self.base_dir / date

        if not date_dir.exists():
            print(f"ℹ No experiments found for {date}")
            return []

        experiment_folders = [folder for folder in date_dir.iterdir() if folder.is_dir()]
        return sorted(experiment_folders)

    def get_next_folder_number(self, sequence: str, date: Optional[str] = None) -> str:
        """
        Get next available folder number for a sequence.

        Args:
            sequence: Sequence type (e.g., "esr", "rabi")
            date: Optional date string. If None, uses today's date.

        Returns:
            Next folder number as string (e.g., "001", "002")

        Example:
            next_num = save_manager.get_next_folder_number("esr")
            # If esr_001 and esr_002 exist, returns "003"
        """
        if date is None:
            date = time.strftime("%Y-%m-%d", time.localtime())

        date_dir = self.base_dir / date

        if not date_dir.exists():
            return "001"

        # Find existing folders for this sequence
        existing = [folder.name for folder in date_dir.iterdir()
                   if folder.is_dir() and folder.name.startswith(f"{sequence}_")]

        if not existing:
            return "001"

        # Extract numbers and find max
        numbers = []
        for folder_name in existing:
            try:
                num = int(folder_name.split('_')[-1])
                numbers.append(num)
            except ValueError:
                continue

        if numbers:
            next_num = max(numbers) + 1
        else:
            next_num = 1

        return f"{next_num:03d}"


if __name__ == "__main__":
    print("\n" + "="*70)
    print(" "*20 + "Save Manager Module Loaded")
    print("="*70)
    print("\nFeatures:")
    print("  - NumPy array storage (.npy)")
    print("  - Separate parameter files (YAML/JSON)")
    print("  - Automatic folder creation: Saved_Data/YYYY-MM-DD/sequence_XXX/")
    print("  - NumPy type conversion for serialization")
    print("\nExample usage:")
    print("  save_manager = SaveManager()")
    print("  save_manager.save_experiment(data, params, '001', 'esr')")
    print("  data, params = save_manager.load_experiment('path/to/data_001.npy')")
    print("="*70 + "\n")
