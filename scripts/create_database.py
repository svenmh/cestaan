import pandas as pd
import sqlite3
import os
import logging
from collections import Counter

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

source_dir = "source_data"

def normalize_headers(headers):
    seen = Counter()
    clean_headers = []
    for col in headers:
        norm = col.strip().replace("/", "_").replace(" ", "_").upper()
        seen[norm] += 1
        if seen[norm] > 1:
            norm = f"{norm}_{seen[norm]}"
        clean_headers.append(norm)
    return clean_headers

def excel_to_sqlite(output_dir: str) -> sqlite3.Connection:
    """
    Converts multiple sheets from an Excel file into separate tables in an SQLite database.

    Parameters:
    - output_dir (str): Directory where the SQLite database will be saved.

    Returns:
    - sqlite3.Connection: An open SQLite connection instance.
    """
    
    def ensure_directory_exists(directory: str) -> None:
        if not os.path.exists(directory):
            os.makedirs(directory)

    def create_database(sheet_names: list, database_name: str, output_dir: str) -> sqlite3.Connection:
        ensure_directory_exists(output_dir)
        database_path = os.path.join(output_dir, database_name)

        # Remove existing database if it exists
        if os.path.exists(database_path):
            os.remove(database_path)
            logging.info(f"Removed existing database: {database_path}")

        # Connect to the database (or create it if it doesn't exist)
        conn = sqlite3.connect(database_path, check_same_thread=False)
        cursor = conn.cursor()
        logging.info(f"Connected to database: {database_path}")

        for sheet_name in sheet_names:
            logging.info(f"Processing sheet: {sheet_name}")
            # Read the Excel file to get the data for each sheet
            df = pd.read_excel(f'{source_dir}/n2_thresholds.xlsx', sheet_name=sheet_name, index_col=0)
            df.index.name = 'Gene'
            df.reset_index(inplace=True)
            headers = df.columns.tolist()
            #Clean up headers
            clean_headers = normalize_headers(headers)
            df.columns = clean_headers
            headers = clean_headers

            # Create the table with dynamic column names
            column_defs = [f'"{header}" REAL' if header.lower() != 'gene' else f'"{header}" TEXT' for header in headers]
            create_table_sql = f"""
            CREATE TABLE IF NOT EXISTS n2_threshold_{sheet_name.split()[-1]} (
                {', '.join(column_defs)}
            )
            """
            cursor.execute(create_table_sql)
            logging.info(f"Created table for n2_threshold_{sheet_name.split()[-1]}")

            # Insert the data
            for row in df.itertuples(index=False, name=None):
                placeholders = ', '.join(['?' for _ in row])
                insert_sql = f"INSERT INTO n2_threshold_{sheet_name.split()[-1]} VALUES ({placeholders})"
                cursor.execute(insert_sql, row)

        # Also add daf2 thresholds
        daf2_sheet_names = ['Threshold 1', 'Threshold 2', 'Threshold 3', 'Threshold 4']
        for sheet_name in daf2_sheet_names:
            logging.info(f"Processing daf2 sheet: {sheet_name}")
            df = pd.read_excel(f'{source_dir}/daf2_thresholds.xlsx', sheet_name=sheet_name, index_col=0)
            df.index.name = 'Gene'
            df.reset_index(inplace=True)
            headers = df.columns.tolist()
            #Clean up headers
            clean_headers = normalize_headers(headers)
            df.columns = clean_headers
            headers = clean_headers

            # Create the table with dynamic column names
            column_defs = [f'"{header}" REAL' if header.lower() != 'gene' else f'"{header}" TEXT' for header in headers]
            create_table_sql = f"""
            CREATE TABLE IF NOT EXISTS daf2_threshold_{sheet_name.split()[-1]} (
                {', '.join(column_defs)}
            )
            """
            cursor.execute(create_table_sql)
            logging.info(f"Created table for daf2_threshold_{sheet_name.split()[-1]}")

            # Insert the data
            for row in df.itertuples(index=False, name=None):
                placeholders = ', '.join(['?' for _ in row])
                insert_sql = f"INSERT INTO daf2_threshold_{sheet_name.split()[-1]} VALUES ({placeholders})"
                cursor.execute(insert_sql, row)

        # Also add Male thresholds
        Male_sheet_names = ['Threshold 1', 'Threshold 2', 'Threshold 3', 'Threshold 4']
        for sheet_name in Male_sheet_names:
            logging.info(f"Processing Male sheet: {sheet_name}")
            df = pd.read_excel(f'{source_dir}/Male_thresholds.xlsx', sheet_name=sheet_name, index_col=0)
            df.index.name = 'Gene'
            df.reset_index(inplace=True)
            headers = df.columns.tolist()
            #Clean up headers
            clean_headers = normalize_headers(headers)

            header_counts = Counter(headers)
            for h, count in header_counts.items():
                if count > 1:
                    print(f"⚠️  Duplicate normalized column: {h} (appears {count} times)")
            df.columns = clean_headers
            headers = clean_headers

            # Create the table with dynamic column names
            column_defs = [f'"{header}" REAL' if header.lower() != 'gene' else f'"{header}" TEXT' for header in headers]            
            create_table_sql = f"""
            CREATE TABLE IF NOT EXISTS Male_threshold_{sheet_name.split()[-1]} (
                {', '.join(column_defs)}
            )
            """
            cursor.execute(create_table_sql)
            logging.info(f"Created table for Male_threshold_{sheet_name.split()[-1]}")

            # Insert the data
            for row in df.itertuples(index=False, name=None):
                placeholders = ', '.join(['?' for _ in row])
                insert_sql = f"INSERT INTO Male_threshold_{sheet_name.split()[-1]} VALUES ({placeholders})"
                cursor.execute(insert_sql, row)

        # Also add Herm thresholds
        Herm_sheet_names = ['Threshold 1', 'Threshold 2', 'Threshold 3', 'Threshold 4']
        for sheet_name in Herm_sheet_names:
            logging.info(f"Processing Herm sheet: {sheet_name}")
            df = pd.read_excel(f'{source_dir}/Herm_thresholds.xlsx', sheet_name=sheet_name, index_col=0)
            df.index.name = 'Gene'
            df.reset_index(inplace=True)
            headers = df.columns.tolist()
            #Clean up headers
            clean_headers = normalize_headers(headers)
            df.columns = clean_headers
            headers = clean_headers

            # Create the table with dynamic column names
            column_defs = [f'"{header}" REAL' if header.lower() != 'gene' else f'"{header}" TEXT' for header in headers]            
            create_table_sql = f"""
            CREATE TABLE IF NOT EXISTS Herm_threshold_{sheet_name.split()[-1]} (
                {', '.join(column_defs)}
            )
            """
            cursor.execute(create_table_sql)
            logging.info(f"Created table for Herm_threshold_{sheet_name.split()[-1]}")

            # Insert the data
            for row in df.itertuples(index=False, name=None):
                placeholders = ', '.join(['?' for _ in row])
                insert_sql = f"INSERT INTO Herm_threshold_{sheet_name.split()[-1]} VALUES ({placeholders})"
                cursor.execute(insert_sql, row)

        # Add additional files
        additional_files = ['daf2_avg.xlsx', 'n2_avg.xlsx', 'daf2_percent.xlsx', 'n2_percent.xlsx','Male_avg.xlsx','Male_pct.xlsx','Herm_avg.xlsx','Herm_pct.xlsx']
        for file_name in additional_files:
            logging.info(f"Processing additional file: {file_name}")
            df = pd.read_excel(f'{source_dir}/{file_name}', index_col=0)
            df.index.name = 'Gene'
            df.reset_index(inplace=True) # Each has one sheet named 'Sheet1'

            df = df.loc[:, ~df.columns.str.lower().isin(['total'])]

            headers = df.columns.tolist()

            # Create the table with dynamic column names
            table_name = file_name.split('.')[0]  # Use the filename without extension as table name
            column_defs = [f'"{header}" REAL' if header.lower() != 'gene' else f'"{header}" TEXT' for header in headers]
            create_table_sql = f"""
            CREATE TABLE IF NOT EXISTS {table_name} (
                {', '.join(column_defs)}
            )
            """
            cursor.execute(create_table_sql)
            logging.info(f"Created table for {table_name}")

            # Insert the data
            for row in df.itertuples(index=False, name=None):
                placeholders = ', '.join(['?' for _ in row])
                insert_sql = f"INSERT INTO {table_name} VALUES ({placeholders})"
                cursor.execute(insert_sql, row)

        # Commit the changes but do not close the connection
        conn.commit()
        logging.info(f"Database '{database_path}' created successfully with tables for sheets: {', '.join(sheet_names)} and daf2 thresholds, plus additional files.")
        return conn

    # Define the sheet names
    sheet_names = ['Threshold 1', 'Threshold 2', 'Threshold 3', 'Threshold 4']

    # Create SQLite database from the sheets and get the connection instance
    conn: sqlite3.Connection = create_database(sheet_names, 'output.sqlite', output_dir)

    # Return the SQLite connection instance
    return conn

if __name__ == "__main__":
    output_dir = "data"  # Hardcoded output directory
    excel_to_sqlite(output_dir)
