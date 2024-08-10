# FMIndexR

This package allows the user to efficiently build an FM index from an input FASTA file. Data structures of FM index can also be saved in a user-defined output folder.

The obtained FM index can then be used to look for patterns inside a sequence, using implementation of Backward search algorithm.

## FM index creation

Using the FM_index_from_FASTA function the FM index 
can easily be created. The function requires different parameters:

-   **input** = the input path to the FASTA file

-   **output** = the output path to the folder in which to save the data structures

-   **save** = regulates whether or not to effectively save the data structures to the output folder (default = TRUE)

```
input_file <- system.file("extdata", "prova.txt", package = "FMIndexR")
output_path <- system.file("output", package = "FMIndexR")
FM_index <- FM_index_from_FASTA(input_file, output_path, save = FALSE)
```

The FM index is a list-like object of class "FM_index" and contains:

-   **Suffix array** : a data structure containing all suffixes of the input string and its associated indexes ordered alphabetically

-   **Burrows-Wheeler Transform (BWT)**: a compression of the original string, generated starting from the suffix array 
for efficiency (compared to the procedure involving a matrix). 
The BWT tends to cluster together similar characters.

The BWT can be built starting from the suffix array, by taking, for each index of the suffix array, the character of the original sequence corresponding to the i-th index - 1

-   **Occurrencies matrix (Occ)**: a table showing the cumulative count of characters of the BWT

-   **Count array ( C )** : for each character in the sorted BWT, the number of characters that are alphabetically smaller than the current character

-   Additionally, also the **FASTA header** is stored 
        to preserve information regarding the original sequence

## Pattern Search

The algorithm I implemented to perform pattern search is 
    Backward Search algorithm, that looks for patterns inside a sequence  using the FM index, by looking at the pattern starting from its end.

Patterns can be searched using the BackwardSearch function, that takes as input:

-   **FM_index** = an object of class FM_index, obtained using the previous function (FM_index_from_FASTA)

-   **pattern** = a string containing a pattern to look for

-   **store_elems** = regulates whether or not to return the original sequence, the index of the pattern occurrencies and the pattern (default = FALSE). This can be done in case a manual check wants to be made.

-   **NOTE**: if no pattern is found, NULL is returned

```{r}
result <- BackwardSearch(FM_index, 'CC', TRUE)
```
The implementation was made following an online explanation that can be found [here](https://tinyurl.com/bwt-reference)

The idea is that the algorithm returns a range of the suffix array in which suffixes starting with the pattern are present. 
The range is initialized with start = 1 (second row of the suffix array, since the first row is always \$) and end = length of the sequence - 1.

At each iteration the following steps are performed:

-   start = Occ[start - 1, letter] + C[letter]
-   end = Occ[end, letter] + C[letter] - 1

If start \> end, it means that the pattern is not found, so the iterations stop. This first implementation could not account for exact matches inside the suffix array 

(i.e those suffixes that entirely matched the pattern) 

In order to fix that, a small modification has been made: the algorithm first looks for exact matches inside the suffix array. 
If exact matches are found, the start of the final range is updated with the value returned by the exact match, if present. The iteration is stopped if start \> stop AND if no exact match is met

>Made -with love- by Gabriele Oliveto
