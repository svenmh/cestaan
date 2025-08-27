import sqlite3
import pandas as pd

db_conn = None

def init():
    global db_conn
    db_conn = sqlite3.connect("data/output.sqlite", check_same_thread=False)

def get_cell_types(genotype: str):
    """
    Retrieve column names from the database, excluding the 'Gene' column.
    """
    cursor = db_conn.cursor()

    # Get the column names
    cursor.execute(f"PRAGMA table_info({genotype}_threshold_1)")
    columns = cursor.fetchall()

    # Exclude the "Gene" column
    column_names = [col[1] for col in columns if col[1].lower() != 'gene' and "total" not in col[1].lower()]

    return column_names

def get_gene_types(genotype: str, threshold: str) -> list[str]:
    cursor = db_conn.cursor()
    
    # Construct the table name based on the threshold
    table_name = f"{genotype}_threshold_{threshold}"

    # Query to retrieve all genes for the specified threshold
    query = f"SELECT Gene FROM {table_name}"
    cursor.execute(query)
    results = cursor.fetchall()

    return [r[0] for r in results]

def get_avg_percent_for_genes(genotype: str, gene_list: list) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Retrieve average expression levels and proportions for a list of genes from the database,
    including support for wildcards (genes ending with a star).

    Parameters:
    - genotype (str): The genotype to use for the table names.
    - gene_list (list): A list of gene names, which may include names ending with a star.

    Returns:
    - tuple: A tuple containing two DataFrames (avg_df, percent_df).
    """
    cursor = db_conn.cursor()
    
    # Prepare the queries for average and percent tables
    avg_query_parts = []
    percent_query_parts = []
    parameters = []
    
    for gene in gene_list:
        if gene.endswith('*'):
            # If the gene ends with a star, use LIKE for wildcard matching
            avg_query_parts.append("Gene LIKE ?")
            percent_query_parts.append("Gene LIKE ?")
            parameters.append(gene[:-1] + '%')  # Replace '*' with '%'
        else:
            avg_query_parts.append("Gene = ?")
            percent_query_parts.append("Gene = ?")
            parameters.append(gene)

    # Combine the query parts with OR
    avg_query = f"SELECT * FROM {genotype}_avg WHERE {' OR '.join(avg_query_parts)}"
    percent_query = f"SELECT * FROM {genotype}_percent WHERE {' OR '.join(percent_query_parts)}"
    
    print(avg_query)
    print(percent_query)

    # Execute the queries
    cursor.execute(avg_query, parameters)
    avg_results = cursor.fetchall()
    
    cursor.execute(percent_query, parameters)
    percent_results = cursor.fetchall()

    print(sorted([r[0] for r in avg_results]))
    print("......")
    print(sorted([r[0] for r in percent_results]))

    # Get cell types for DataFrame columns
    cell_types = get_cell_types(genotype)

    # Create DataFrames for average and percent
    avg_df = pd.DataFrame(avg_results, columns=['Gene', *cell_types])
    percent_df = pd.DataFrame(percent_results, columns=['Gene', *cell_types])

    return avg_df, percent_df  # Return as a tuple of DataFrames

def get_cell_type_by_gene(genotype: str, threshold: str, gene: str = None):
    """
    Retrieve cell types and their expression levels for a specific threshold and gene from the database.
    
    Parameters:
    - threshold (str): The threshold to use for the table name.
    - gene (str, optional): The gene to look up. If provided, filters results by this gene.
    
    Returns:
    - pd.DataFrame: A DataFrame with cell types as rows and their expression levels for the specified gene.
    """
    cursor = db_conn.cursor()
    
    # Construct the table name based on the threshold
    table_name = f"{genotype}_threshold_{threshold}"

    if gene:
        # If a gene is provided, filter the results based on the gene
        query = f"SELECT * FROM {table_name} WHERE Gene = ?"
        cursor.execute(query, (gene,))
        results = cursor.fetchall()

        if len(results) == 0:
            return None

        results = results[0][1:]

        # Create a DataFrame with cell types as rows and expression levels as values
        column_names = get_cell_types(genotype)
        # Create a DataFrame with cell types as rows
        df = pd.DataFrame(sorted(list(zip(column_names, results)), key=lambda x: x[1], reverse=True), index=column_names, columns=['Cell Type', 'Expression Level'])
        df.index.name = 'Cell Type'
        
        return df  # Return as DataFrame
    
    return pd.DataFrame()  # Return an empty DataFrame if no gene is provided

def get_gene_by_cell_type(genotype: str, threshold: str, cell_type: str):
    """
    Retrieve genes and their expression levels for a specific cell type from the database,
    sorted from highest to lowest expression level.
    
    Parameters:
    - threshold (str): The threshold to use for the table name.
    - cell_type (str): The cell type to look up.
    
    Returns:
    - pd.DataFrame: A DataFrame with genes as rows and their expression levels for the specified cell type.
    """
    import pandas as pd
    cursor = db_conn.cursor()
    
    # Construct the table name based on the threshold
    table_name = f"{genotype}_threshold_{threshold}"

    # Query to get genes and their expression levels for the specified cell type
    query = f"SELECT Gene, \"{cell_type}\" FROM {table_name} WHERE \"{cell_type}\" IS NOT NULL ORDER BY \"{cell_type}\" DESC"
    cursor.execute(query)
    results = cursor.fetchall()

    # Create a DataFrame with genes as rows and their expression levels
    df = pd.DataFrame(results, columns=['Gene', 'Expression Level'])
    
    return df  # Return as DataFrame

def get_genes_bulk(genotype: str, threshold: str, gene_list: list) -> pd.DataFrame:
    """
    Retrieve genes and their expression levels for a list of genes from the database,
    including support for wildcards (genes ending with a star).

    Parameters:
    - genotype (str): The genotype to use for the table name.
    - threshold (str): The threshold to use for the table name.
    - gene_list (list): A list of gene names, which may include names ending with a star.

    Returns:
    - pd.DataFrame: A DataFrame with genes as rows and their expression levels.
    """
    import pandas as pd
    cursor = db_conn.cursor()
    
    # Construct the table name based on the threshold
    table_name = f"{genotype}_threshold_{threshold}"

    # Prepare the query with placeholders
    query_parts = []
    parameters = []
    
    for gene in gene_list:
        if gene.endswith('*'):
            # If the gene ends with a star, use LIKE for wildcard matching
            query_parts.append("Gene LIKE ?")
            parameters.append(gene[:-1] + '%')  # Replace '*' with '%'
        else:
            query_parts.append("Gene = ?")
            parameters.append(gene)

    # Combine the query parts with OR
    query = f"SELECT * FROM {table_name} WHERE {' OR '.join(query_parts)}"
    print(query)
    cursor.execute(query, parameters)
    results = cursor.fetchall()

    print(results)

    cell_types = get_cell_types(genotype)

    # Create a DataFrame with genes as rows and their expression levels
    df = pd.DataFrame(results, columns=['Gene', *cell_types])
    
    return df  # Return as DataFrame

def df_to_csv(df: pd.DataFrame) -> str:
    """
    Convert a DataFrame to a CSV string for download.

    Parameters:
    - df (pd.DataFrame): The DataFrame to convert.

    Returns:
    - str: The CSV string representation of the DataFrame.
    """
    return df.to_csv(index=False)
