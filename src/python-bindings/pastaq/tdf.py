import sqlite3
from typing import List, Dict
import graphviz

class tdf:
    def __init__(self, db_path: str):
        """
        Initialize the SQLiteExplorer with the path to the SQLite database.
        :param db_path: Path to the SQLite database file.
        """
        self.db_path = db_path
        self.connection = None

    def connect(self):
        """Establish a connection to the SQLite database."""
        self.connection = sqlite3.connect(self.db_path)

    def close(self):
        """Close the connection to the SQLite database."""
        if self.connection:
            self.connection.close()
            self.connection = None

    def get_tables(self) -> List[str]:
        """
        Retrieve the list of tables in the database.
        :return: List of table names.
        """
        self._ensure_connection()
        cursor = self.connection.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = [row[0] for row in cursor.fetchall()]
        return tables

    def get_schema(self, table_name: str) -> None:
        """
        Retrieve and pretty-print the schema of a specified table.
    
        :param table_name: Name of the table.
        """
        self._ensure_connection()
        cursor = self.connection.cursor()
    
        # Fetch the schema information
        cursor.execute(f"PRAGMA table_info({table_name});")
        schema = cursor.fetchall()
    
        if not schema:
            print(f"Table '{table_name}' does not exist or has no schema.")
            return

        # Pretty-print the schema
        print(f"\nSchema for table '{table_name}':")
        print(f"{'Column ID':<10} {'Column Name':<20} {'Data Type':<15} {'Not Null':<10} {'Default Value':<15} {'Primary Key':<10}")
        print("-" * 80)

        for column in schema:
            column_id, name, data_type, not_null, default_value, primary_key = column
            print(f"{column_id:<10} {name:<20} {data_type:<15} {not_null:<10} {str(default_value):<15} {primary_key:<10}")
        print("-" * 80)

    def query(self, sql_query: str) -> List[Dict]:
        """
        Execute an SQL query and return the results as a list of dictionaries.
        :param sql_query: The SQL query string.
        :return: Query results as a list of dictionaries.
        """
        self._ensure_connection()
        cursor = self.connection.cursor()
        cursor.execute(sql_query)
        columns = [desc[0] for desc in cursor.description]
        results = [dict(zip(columns, row)) for row in cursor.fetchall()]
        return results

    def relationships(self):
        """
        Print the relationships (foreign keys) between tables in the SQLite database.
        """
        self._ensure_connection()
        cursor = self.connection.cursor()

        # Retrieve the list of tables in the database
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        tables = cursor.fetchall()

        # Initialize a dictionary to store relationships
        relationships = {}

        # Iterate through all tables to identify relationships
        for table in tables:
            table_name = table[0]
            cursor.execute(f"PRAGMA foreign_key_list('{table_name}');")
            fk_info = cursor.fetchall()

            # Store relationships in the dictionary
            for fk in fk_info:
                source_table = fk[3]  # Referenced table
                target_table = table_name  # Table containing the foreign key
                column_relationship = f"{fk[4]} -> {fk[3]}"  # Foreign key column -> Referenced column

                if source_table not in relationships:
                    relationships[source_table] = []
                relationships[source_table].append((target_table, column_relationship))

        # Print relationships in a readable format
        print("Table Relationships (Foreign Keys):\n")
        for source_table, targets in relationships.items():
            print(f"Table '{source_table}' is referenced by:")
            for target_table, column_relationship in targets:
                print(f"  -> {target_table} ({column_relationship})")
            print()


    def _ensure_connection(self):
        """Ensure that the database connection is established."""
        if not self.connection:
            raise RuntimeError("Database connection is not established. Call `connect()` first.")
