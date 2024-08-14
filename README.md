# FMIndexR

This package allows the user to efficiently build an FM index from an input DNA FASTA file. Data structures composing the FM index can also be saved as stand-alone files in a user-defined output folder.

The obtained FM index can then be used to look for patterns inside a sequence, using an implementation of Backward search algorithm.

## **_<u>FM index creation  </u>_**

Using the **_FM_index_from_FASTA_** function an FM index 
can be easily created. The function requires different parameters:

-   **_input =_**  the input path to the FASTA file

-   **_output =_**  the output path to the folder in which to save the data structures

-   **_save =_**  regulates whether or not to effectively save the data structures to the output folder (**default = TRUE**)

```
input_file <- system.file("extdata", "prova_vignette.txt", package = "FMIndexR")
output_path <- system.file("output", package = "FMIndexR")
FM_index <- FM_index_from_FASTA(input_file, output_path, save = FALSE)
```

The FM index is a list-like object of class "FM_index" and contains:

-   **_Suffix array :_**  a data structure containing all suffixes of the input string and its associated starting indexes, ordered alphabetically

-   **_Burrows-Wheeler Transform (BWT) :_** a compression of the original string, generated starting from the suffix array for efficiency (compared to the procedure involving a matrix). The BWT tends to cluster together similar characters. 

**_NOTE:_** The BWT can be built starting from the suffix array, by taking, for each index i of the suffix array, the character of the original sequence corresponding to the i-th index - 1

-   **_Occurrencies matrix (Occ) :_** a table showing the cumulative count of characters of the BWT

-   **_Count array ( C ) :_**  for each character in the sorted BWT, the number of characters that are alphabetically smaller than the current character. The count array can be viewed also as the first index of occurrence of each character inside the sorted BWT

-   Additionally, also the **_FASTA header_** is stored to preserve information regarding the original sequence

## **_<u>Pattern Search</u>_**

The algorithm I used to perform pattern search is an implementation of the Backward Search algorithm, that looks for patterns inside a sequence using the FM index, by looking at the pattern starting from its end.

Patterns can be searched using the **_BackwardSearch_** function, that takes as input:

-   **_FM_index =_**  an object of class FM_index, obtained using the previous function (**_FM_index_from_FASTA_**)

-   **_pattern =_**  a non-empty string containing a pattern to look for

-   **_store_elems =_**  regulates whether or not to return the original sequence, the index of the pattern occurrencies and the pattern (**default = FALSE**). This can be done in case a manual check wants to be made.

The function shows the user the number of found patterns and their indexes. **If no pattern is found, NULL is returned**

```{r}
result <- BackwardSearch(FM_index, 'GATG', TRUE)
```
The implementation was built following an online explanation that can be found [here](http://blog.thegrandlocus.com/2016/07/a-tutorial-on-burrows-wheeler-indexing-methods)

The main idea is that the algorithm returns a range of rows in the suffix array in which suffixes starting with the pattern are present.
    
Occurrencies are then found inside the original sequence by looking at the indexes associated with the suffix inside the suffix array

The range is initialized with:
- **_start =_**  1 (second row of the suffix array, since the first row is always \$)
- **_end =_**  length of the BWT - 1.

At each iteration (going backwards along the pattern) the following steps are performed:
    
-   for the first iteration only **_start =_** C[letter]

-   then **_start =_**  Occ[start - 1, letter] + C[letter]

-   **_end =_**  Occ[end, letter] + C[letter] - 1

If **_start \> end_**, it means that the **_pattern is not found_**, so the iterations stops. 

## _**<u>Additional Notes</u>**_

- **R version :**  4.3.2
- **Dependencies :**
	- **Biostrings :**  2.70.3
	- **IRanges :**  2.36.0
	- **zoo :**  1.8-12

>Made -with love- by Gabriele Oliveto
