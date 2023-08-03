from pandas import DataFrame
from os.path import join as osp_join
from time import strftime

def export_to_csv(to_export: dict, file_name: str = "", directory: str = "results") -> None:
    """
    Export data from a dictionary to a CSV file.

    Params:
        to_export (dict): A dictionary where keys are column names and values are lists of column values.
                          ALL VALUE LISTS MUST BE OF THE SAME LENGTH TO CREATE A DATAFRAME.
        file_name (str): Optional custom string appended before the timestamp in the resulting filename.
        directory (str): Optional custom target directory. The default value is "results" directory.

    """
    print("Exporting results.")

    df = DataFrame(to_export)
    # Create a CSV file in the specified directory with the given file name and timestamp
    df.to_csv(osp_join(directory, file_name + strftime("%Y_%m_%d-%H_%M_%S") + ".csv")) 

def main():
    # Example use:
    export_to_csv(
        {
            'test_column_1': ['val_1', 'val_2'],
            'test_column_2': ['val_1']
        },
        file_name="test"
    )

if __name__ == '__main__':
    main()
