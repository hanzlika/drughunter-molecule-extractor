from DECIMER import predict_SMILES
import os

def recognize_segments(input_directory, output_directory):
    """Iterates through the input_directory
    and attempts to recognize every file within it
    as a SMILES string, which is then appended onto
    the output recognized_smiles.txt
    """
    with open (f"{output_directory}/recognized_smiles.txt", "a") as file_object:
        for filename in os.listdir(input_directory):
            image_path = os.path.join(input_directory, filename)
          
            SMILES = predict_SMILES(image_path)
            file_object.write(f"{SMILES}\n")

def main():
    input_directory = "segmentation/segments"
    output_directory = "recognition"
    recognize_segments(input_directory, output_directory)

if __name__ == '__main__':
    main()