import sys
import pandas as pd
from sklearn.preprocessing import StandardScaler
import joblib

# Load the normalization and trained model files using joblib
def load_model_and_scaler(scaler_filepath, model_filepath):
    # Load the scaler
    scaler = joblib.load(scaler_filepath)
    # Load the trained model
    model = joblib.load(model_filepath)
    return scaler, model

# Function to load all features from the csv file
def load_features_from_csv(filepath):
    # Load the CSV file, specifying that the first column is an index (if needed)
    data = pd.read_csv(filepath, index_col=0)  # Assuming 'Unnamed: 0' is the index, we ignore it
    return data

# Function to map numeric predictions to labels
def map_predictions(predictions):
    # Mapping numeric predictions to classes
    label_mapping = {0: 'Cardiotoxic', 1: 'Enterotoxic', 2: 'Neurotoxic', 3: 'Non-toxic'}
    mapped_predictions = [label_mapping[pred] for pred in predictions]
    return mapped_predictions

# Predict function
def predict_with_model(scaler_file, model_file, csv_file):
    # Load the normalization scaler and the trained model
    scaler, model = load_model_and_scaler(scaler_file, model_file)

    # Load the features (entire CSV file as features)
    feature_data = load_features_from_csv(csv_file)

    # Normalize the feature data using the scaler
    normalized_features = scaler.transform(feature_data)

    # Use the trained model to make predictions
    predictions = model.predict(normalized_features)

    # Map predictions to class labels
    mapped_predictions = map_predictions(predictions)

    return mapped_predictions

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <scaler.pkl> <model.pkl> <abc.csv>")
        sys.exit(1)

    # Arguments from command line
    scaler_pkl = sys.argv[1]  # Path to your normalization .pkl file (scaler)
    model_pkl = sys.argv[2]  # Path to your trained model .pkl file
    csv_file = sys.argv[3]  # Path to your csv file (abc.csv)

    # Run the prediction
    predictions = predict_with_model(scaler_pkl, model_pkl, csv_file)

    # Create a DataFrame for the results
    output_df = pd.DataFrame(predictions, columns=['Predicted Class'])

    # Save the predictions to a CSV file
    output_df.to_csv('predictions_output.csv', index=False)

    # Output confirmation
    print("Predictions saved to 'predictions_output.csv'")
