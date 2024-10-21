# Import necessary libraries
import argparse
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
import seaborn as sns

def load_data(train_file, test_file):
    """Load train and test datasets from CSV files."""
    train_data = pd.read_csv(train_file)
    test_data = pd.read_csv(test_file)
    
    X_train = train_data.iloc[:, :-1]
    y_train = train_data.iloc[:, -1]
    
    X_test = test_data.iloc[:, :-1]
    y_test = test_data.iloc[:, -1]
    
    return X_train, y_train, X_test, y_test

def train_and_evaluate(X_train, y_train, X_test, y_test):
    """Train a Random Forest model, perform GridSearchCV, and evaluate the model."""
    
    # Initialize Random Forest
    rf = RandomForestClassifier()

    # Define hyperparameter grid
    param_grid = {
        'n_estimators': [50, 100, 200],
        'max_depth': [None, 10, 20, 30],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4]
    }

    # Perform grid search with cross-validation
    grid_search = GridSearchCV(estimator=rf, param_grid=param_grid, cv=10, n_jobs=-1, verbose=2)
    grid_search.fit(X_train, y_train)

    # Get the best model
    best_knn = grid_search.best_estimator_
    print("Best Parameters:", grid_search.best_params_)

    # Make predictions
    y_pred_train = best_knn.predict(X_train)
    y_pred_test = best_knn.predict(X_test)

    # Evaluate the model
    print("Classification Report Training Set:")
    print(classification_report(y_train, y_pred_train))

    print("Classification Report of Test Set:")
    print(classification_report(y_test, y_pred_test))

    # Generate the confusion matrix
    conf_matrix_train = confusion_matrix(y_train, y_pred_train)
    print("Confusion Matrix of Training Set:")
    print(conf_matrix_train)

    conf_matrix_test = confusion_matrix(y_test, y_pred_test)
    print("Confusion Matrix of Test Set:")
    print(conf_matrix_test)

def main():
    """Main function to parse arguments and run the model."""
    parser = argparse.ArgumentParser(description="Train and evaluate a Random Forest classifier using GridSearchCV.")
    
    # Add arguments for train and test CSV files
    parser.add_argument('--train', type=str, required=True, help='Path to the training data CSV file')
    parser.add_argument('--test', type=str, required=True, help='Path to the testing data CSV file')

    # Parse the arguments
    args = parser.parse_args()
    
    # Load the data
    X_train, y_train, X_test, y_test = load_data(args.train, args.test)
    
    # Train the model and evaluate
    train_and_evaluate(X_train, y_train, X_test, y_test)

if __name__ == "__main__":
    main()