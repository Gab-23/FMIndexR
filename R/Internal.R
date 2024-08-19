
#BackwardSearch

.getRangeandCount <- function(rev_pattern_array,Occ,C,BWT){

    getRange <- function(rev_pattern_array,Occ,C,BWT){
        start <- 1
        end <- (nchar(BWT) - 1)
        for (idx in seq_along(rev_pattern_array)) {
            char <- rev_pattern_array[idx]
            if (idx == 1) {
                start <- C[[char]]
                } else {
                    start <- C[[char]] +
                    Occ[as.character(start - 1), char]}
                end <- C[[char]] + Occ[as.character(end), char] - 1
                if (start > end) {
                    return(NULL)
                    break}}

        return(list(start = start, end = end))}

    range <- getRange(rev_pattern_array,Occ,C,BWT)
    if (is.null(range)) {
        return(list(num_pattern = 0))
        } else {
            start <- range$start
            end <- range$end

            final_range <- as.character((start):(end))
            num_pattern <- length(final_range)
            return(list(range = final_range, num_pattern = num_pattern))}}

.getIndexes <- function(BWT,SA,C,Occ,final_range,num_pattern){

    rebuild_idx <- function(BWT,SA,C,Occ,val){
        i <- 0
        while (TRUE) {
            char <- substr(BWT,val+1,val+1)
            val <- C[[char]] + Occ[as.character(val),char] - 1
            i <- i + 1
            if (val %% 2 == 0) {
                idx <- SA[as.character(val),"idx"] + i
                if (idx >= nchar(BWT)) {
                    idx <- idx %% nchar(BWT)}
                return(idx)}}}

    indexes <- vector("list",length(num_pattern))

    for (elem in seq_along(final_range)) {
        num <- final_range[elem]
        if (num %in% rownames(SA)) {
            indexes[elem] <- SA[num,"idx"]
            } else {
                indexes[elem] <- rebuild_idx(BWT,SA,C,Occ,as.integer(num))}}
    return(indexes)}

.reverseBWT <- function(BWT,C,Occ){
    BWT_arr <- strsplit(BWT,split = '')[[1]]
    char<- "$"
    val <- which(BWT_arr == char) - 1
    seq <- character(length(BWT_arr))
    for (i in seq_along(BWT_arr)) {
        seq[i] <- char
        val <- C[[char]] + Occ[as.character(val),char] - 1
        char <- BWT_arr[val+1]}
    seq <- paste(rev(seq), collapse = '')
    seq <- substr(seq,1,nchar(seq)-1)
    return(seq)}

.get_original_sequence <- function(SA,BWT,Occ,C){
    if (0 %in% SA$idx) {
        original_sequence <- SA[SA$idx == 0,"suffix"]
        original_sequence <- substr(original_sequence,
                                    1,nchar(original_sequence)-1)
        } else {
            original_sequence <- .reverseBWT(BWT,C,Occ)}
    return(original_sequence)}

#FM_index_from_FASTA

.save <- function(seq,seq_path,SA,SA_path,BWT,BWT_path,Occ,Occ_path,C,C_path){

    C <- as.data.frame(C)
    writeLines(seq,con = seq_path)
    utils::write.table(SA,file = SA_path,sep = "\t",row.names = FALSE)
    writeLines(BWT,con = BWT_path)
    utils::write.table(Occ,file = Occ_path,sep = "\t",row.names = FALSE)
    utils::write.table(C,file = C_path,sep = "\t",row.names = TRUE)}

.SuffixArray <- function(input_string) {
    if (nchar(input_string) == 0) {
        stop("Empty sequence detected!")
        } else {
            special_char <- "$"
            complete_string <- paste(input_string,
                                        special_char,
                                        sep = "",
                                        collapse = "")

            complete_string_length <- nchar(complete_string)
            idx_array <- 0:(complete_string_length - 1)

            suffix_array <- vapply(X = idx_array, function(i) {
                actual_idx <- i
                used_idx <- i + 1
                suffix <- substr(complete_string,
                                    used_idx,
                                    complete_string_length)}, character(1))

            suffix_df <- data.frame(idx = idx_array,
                                    suffix = suffix_array)

            sorted_suffix_df <- suffix_df[order(suffix_array), ]
            rownames(sorted_suffix_df) <- 0:(nrow(sorted_suffix_df) - 1)
            return(sorted_suffix_df)}
}

.BWTransform <- function(suffix_array) {

    original_sequence <- suffix_array[suffix_array$idx == 0, 2]
    actual_idx <- suffix_array$idx
    used_idx <- ifelse(actual_idx == 0,
                        nchar(original_sequence),
                        actual_idx)

    BWT <- vapply(X = used_idx, function(i) {
        idx <- i
        substr(original_sequence, idx, idx)}, character(1))

    BWT <- paste(BWT, sep = "", collapse = "")
    return(BWT)}

.OccMatrix <- function(bwt) {

    length_bwt <- nchar(bwt)
    bwt_letters <- strsplit(bwt, split = "")[[1]]
    unique_bwt_letters <- unique(bwt_letters)
    sorted_unique_bwt_letters <- sort(unique_bwt_letters)

    Occ_mat <- matrix(data = NA,
                        nrow = length_bwt,
                        ncol = length(sorted_unique_bwt_letters))

    colnames(Occ_mat) <- sorted_unique_bwt_letters

    for (char in colnames(Occ_mat)) {

        idx_vector <- which(bwt_letters == char)
        Occ_mat[idx_vector, char] <- seq_along(idx_vector)

        count_vec <- Occ_mat[, char]
        count_vec <- zoo::na.locf(count_vec, na.rm = FALSE)
        count_vec[is.na(count_vec)] <- 0

        Occ_mat[, char] <- count_vec}

    rownames(Occ_mat) <- 0:(nrow(Occ_mat) - 1)
    return(Occ_mat)}

.CountArray <- function(bwt) {

    bwt_char <- sort(unique((strsplit(bwt, split = "")[[1]])))
    sorted_bwt <- sort(strsplit(bwt, split = "")[[1]])
    count_vec <- numeric(length(bwt_char))
    names(count_vec) <- bwt_char

    for (char in bwt_char) {
        char_idx <- which(sorted_bwt == char)[1] - 1
        count_vec[char] <- char_idx}

    return(count_vec)}

.FMIndex <- function(sequence_name,SA,BWT,Occ,C){
    FM_index <- list(SequenceName = sequence_name,
                    SuffixArray = SA,
                    BWT = BWT,
                    Occ = Occ,
                    CountArray = C)

    class(FM_index) <- "FM_index"
    return(FM_index)
}



