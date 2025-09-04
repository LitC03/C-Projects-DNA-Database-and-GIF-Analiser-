# DNA Database Manager

The DNA Database Manager is a command-line utility, developed in C++ (C++11 standard), and engineered for the management and manipulation of DNA sequences derived from FASTA (`.fna`) formatted files. This application facilitates the loading of multiple DNA sequences into an in-memory database, enabling a comprehensive suite of search and modification operations, with the capability to persist the results. The program leverages a custom, linked-list-based data structure to ensure efficient processing and handling of substantial DNA sequence data.

## Core Functionality

* **Sequence Loading:** Ingests one or more DNA sequences from `.fna` files into the database.

* **Sequence Processing:** Provides a suite of operations for a selected DNA sequence:

  * **Subsequence Search:** Locates a DNA subsequence using either direct textual input or an external FASTA file as the query.

  * **Subsequence Insertion:** Inserts a new DNA subsequence at a user-specified position from either textual input or a file.

  * **Subsequence Deletion:** Excises a segment of the DNA sequence based on a specified starting position and length.

  * **Subsequence Replacement:** Overwrites a segment of the DNA sequence with a new subsequence provided via textual input or a file.

  * **Result Persistence:** Saves a modified DNA sequence to a new `.fna` file.

* **Database-Wide Analysis:** Conducts a comprehensive search for a specific DNA sequence, sourced from a file, across all sequences currently resident in the database.

* **Performance Measurement:** Quantifies and reports the execution time for file loading and data processing tasks.

## Compilation and Execution

### System Requirements

* A standard C++ compiler (e.g., MinGW).

* The `tictoc.h` header file must be located within the same directory as the source files.

### Compilation & Execution (MinGW)

To compile the application, navigate to the source directory via a command-line interface and execute the following command:
```
g++ DNA_Database_Manager.cpp -o DNA_Database_Manager –std=c++11 
```
To run the compiled program (already in the repo), use the following command:
```
DNA_Database_Manager.exe
```
## Operational Workflow

Upon execution, the application presents an interactive main menu with the following options:

1. **Load DNA(s):**

   * Prompts the user to input the filenames of the `.fna` files for loading. Multiple filenames should be separated by a comma (e.g., `file1.fna, file2.fna`).

2. **Process a DNA:**

   * Displays a list of the DNA files currently loaded in the database.

   * The user selects a file for processing by entering its corresponding numerical index.

   * A submenu is then presented, offering operations such as finding, adding, deleting, replacing, or saving the selected DNA sequence.

3. **Analyse the DNA database:**

   * Prompts the user for the name of an `.fna` file.

   * The system then searches for the sequence contained within that file across all DNA sequences in the database, reporting any and all matches found.

4. **Quit:**

   * Terminates the application.

## Architectural Design

The application's architecture is founded upon a hierarchical data structure implemented through a set of custom classes, forming a nested linked-list configuration:

* **`DNADatabase`**: This primary class serves as the container for the entire database, functioning as a linked list in which each node is an `Index` object.

* **`Index`**: This class encapsulates a single, complete DNA sequence loaded from a file. It functions as a node within the `DNADatabase` and serves as the head of a linked list composed of `Nucleotide` objects.

* **`Nucleotide`**: This class represents an individual DNA base (e.g., A, C, G, T). These objects are linked sequentially to constitute the complete DNA sequence managed by an `Index` object.

This design facilitates dynamic memory allocation and enables highly efficient insertion and deletion operations, thereby avoiding the performance overhead associated with reallocating large, contiguous memory blocks.
