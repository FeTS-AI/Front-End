from abc import ABC, abstractmethod
import pandas as pd


class DatasetStage(ABC):
    @abstractmethod
    def should_run(self, report: pd.DataFrame) -> bool:
        """Establishes if this step should be executed

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): Dataframe containing the current state of the preparation flow

        Returns:
            bool: wether this stage should be executed
        """

    @abstractmethod
    def execute(self, report: pd.DataFrame) -> pd.DataFrame:
        """Executes the stage

        Args:
            index (Union[str, int]): case index in the report
            report (pd.DataFrame): DataFrame containing the current state of the preparation flow

        Returns:
            pd.DataFrame: Updated report dataframe
        """
