# FMIndexR

This package allows the user to efficiently build an FM index from an input DNA FASTA file. Data structures composing the FM index can also 
be saved as stand-alone files in a user-defined output folder.

The obtained FM index can then be used to look for patterns inside a sequence, using an implementation of Backward search algorithm.

**NOTE:** since the package deals with (possibly) large tabular formats, in order to avoid indexing errors the scientific notation has been disabled:
- Both functions will initially set the option scipen = 999
- The standard setting using scipen = 0 will be reset after each usage

## <u>FM index creation </u>

Using the ***FM_index_from_FASTA*** function  an FM index can be easily created. The function requires different parameters:

-   ***input =*** the input path to the FASTA file

-   ***output =*** the output path to the folder in which to save the data structures

-   ***save =*** regulates whether or not to effectively save the data structures to the output folder (**default = TRUE**)

-   ***compress =*** regulates whether or not to compress the FM index by downsampling the suffix array (**default = FALSE**)

```{r}
FM_index <- FM_index_from_FASTA(input_file, 
                                output_path, 
                                save = FALSE, 
                                compress = FALSE)
```
The FM index is a list-like object of class "FM_index" and contains:

-   ***Suffix array :*** a dataframe containing all (or a part of the) 
suffixes of the input string and its associated starting indexes, 
ordered alphabetically.

**NOTE:**  when ***compress = TRUE*** , the suffix array gets compressed, 
retaining one suffix every 32. As the suffix array gets downsampled,
the query time becomes longer. 

In order to understand the optimal compression rate,  a plot has been shown. 

![Image currently not available :(](inst/extdata/size-time-plot.png)

- On the x-axis you find the size in bytes of the suffix array.
- On the y-axis you find the average (n = 30) query time in seconds.

As it appears, keeping one suffix every 32 is the compression rate that grants the best trade-off between query time and suffix array size reduction.

-   ***Burrows-Wheeler Transform (BWT) :*** a compression of 
the original string, that tends to cluster together similar characters. 
It is generated starting from the suffix array for efficiency (compared to the procedure involving a matrix).

    ***NOTE:*** The BWT can be built starting from the suffix array, 
by taking, for each index i of the suffix array,  the character of the original sequence corresponding to the i-th index - 1.

-   ***Occurrencies matrix (Occ) :*** a matrix showing the cumulative count of characters of the BWT.

-   ***Count array ( C ) :*** a vector showing, for each character in the sorted BWT, the number of characters that are alphabetically smaller than the current character. The count array can be viewed also as the first index of occurrence of each character in the sorted BWT

-   Additionally, also the ***FASTA header*** is stored to preserve information regarding the original sequence.

## <u>Pattern Search</u>

The algorithm I used to perform pattern search is an implementation 
of the Backward Search algorithm, that looks for patterns 
inside a sequence using the FM index, 
by looking at the pattern starting from its end.

Patterns can be searched using the ***BackwardSearch*** function, 
that takes as input:

-   ***FM_index =*** an object of class FM_index, 
    obtained using the function ***FM_index_from_FASTA***

-   ***pattern =*** a non-empty string containing the pattern to look for

-   ***store_elems =*** regulates whether or not to return a more detailed 
    output containing the original sequence, the index of the pattern 
    occurrencies and the pattern (**default = FALSE**). 
    This can be done in case a manual check wants to be made

The function shows the user the number of found patterns and their indexes  (which are also returned). **If no pattern is found, NULL is returned**.

```{r}
result <- BackwardSearch(FM_index, 'GATG', TRUE)
```
**NOTE:** Pattern search performed on the uncompressed FM index retrieve the same results as the compressed one.

The implementation was built following an online explanation that can be found [here](http://blog.thegrandlocus.com/2016/07/a-tutorial-on-burrows-wheeler-indexing-methods).

The main idea is that the algorithm returns a range of rows  in the suffix array in which suffixes starting with the pattern are present.

The range is initialized with:

-   ***start =*** 1 (second row of the suffix array,  since the first row always contains \$)
-   ***end =*** length of the BWT - 1.

At each iteration (going backwards along the pattern) the following steps are performed:

-   for the first iteration only ***start =*** C[letter]
-   then ***start =*** Occ[start - 1, letter] + C[letter]

-   ***end =*** Occ[end, letter] + C[letter] - 1

If ***start \> end***, it means that the ***pattern is not found***,  so the iterations stop.

Occurrencies are then found inside the original sequence by looking at the rows pointed by the range. The rows contain the suffix starting with the pattern and the associated index.

The procedure gets more time-consuming and complicated if the suffix array is compressed, since not all elements of the range point 
towards an existing row of the suffix array.

In order to retrieve the index of a missing row a procedure known as **Last to First (LF) mapping** is applied.

**NOTE:** Common implementations of LF mapping use recursion. 
Here, an iterative process has been implemented, hence, some differences between this implementation and the known ones may be present (more on this later).

Given a range R, if the i-th element of the range (**val**)  is not included in the suffix array:

-   The i-th character (**char**) of the BWT is retrieved

-   **val** is updated as follows: C[char] + Occ[char,val] - 1

-   If **val** is a multiple of the compression rate (32) 
then the row corresponding to **val** exists in the suffix array  and the corresponding index in the original string is computed as:

The index pointed by **val** in the suffix array + the number of iterations performed

The main difference with recursive implementations regards out-of-bounds indexes. If the value of the index exceeds the length  of the sequence, then the **modulus operator %%** (index %% length) is used obtain a compatible index.

The difference arises when the **$** symbol is encountered:

- The recursive implementation would return 0 + the number of iterations.

- Here the function returns instead the index associated to $ 
    (always the number of characters - 1) + the number of iterations
    
So for example if **on the first iteration** we encountered the $:

- In the common implementation we would simply return 0
- Here we return the index associated to $ + 1 (iteration) = number of characters - 1 + 1 = number of characters that will be later trimmed to 0 using the %% operator.

Lastly, if **store_elems = TRUE** a detailed output containing the original sequence, the pattern and the indexes is returned. The program will try to look if the original sequence is already present in the suffix array. If not present, the original sequence will be reconstructed from the BWT and that operation may take a longer time.

> Made -with love- by Gabriele Oliveto

