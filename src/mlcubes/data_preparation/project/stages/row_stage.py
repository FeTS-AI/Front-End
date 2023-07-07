from abc import ABC, abstractmethod
from typing import Union
import pandas as pd


class RowStage(ABC):
    @abstractmethod
    def should_run(self, index: Union[str, int], report: pd.DataFrame) -> bool:
        """Establishes if this step should be executed for the given case

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): Dataframe containing the current state of the preparation flow

        Returns:
            bool: wether this stage should be executed
        """

    @abstractmethod
    def execute(self, index: Union[str, int], report: pd.DataFrame) -> pd.DataFrame:
        """Executes the stage on the given case

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
        """
