from abc import ABC, abstractmethod
import pandas as pd
from typing import Tuple


class DatasetStage(ABC):
    @abstractmethod
    def could_run(self, report: pd.DataFrame) -> bool:
        """Establishes if this step could be executed

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): Dataframe containing the current state of the preparation flow

        Returns:
            bool: wether this stage could be executed
        """

    @abstractmethod
    def execute(self, report: pd.DataFrame) -> Tuple[pd.DataFrame, bool]:
        """Executes the stage

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
            bool: Success status
        """

    @abstractmethod
    def get_name(self) -> str:
        """Returns a human readable short name for what the stage does
        Used for printing to the user current status

        Returns:
            str: Stage name
        """
