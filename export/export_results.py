from pandas import DataFrame
from os.path import join as osp_join
from time import strftime

def export_to_csv(to_export:dict, file_name:str = "", directory:str = "results"):
    """
    Params:
        to_export -> dictionary, keys = names of column, values = list of values in columns (ALL VALUE LISTS MUST BE OF THE SAME LENGTH TO CREATE A DATAFRAME)
        file_name -> optional custom string appended before the timestamp in resulting filename
        directory -> optional custom target directory, default value = "results" directory

    """
    print("Exporting results.")

    df = DataFrame(to_export)
    df.to_csv(osp_join(directory, file_name + strftime("%Y_%m_%d-%H_%M_%S") + ".csv")) 

def main():
    export_to_csv(
        dict(
                {
                'test_column_1':
                    ['val_1', 'val_2'],
                'test_column_2':
                    ['val_1']
                }
            ),
        "test"
        )

if __name__ == '__main__':
    main()