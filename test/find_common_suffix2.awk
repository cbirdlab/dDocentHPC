# define the input field separator as a newline character
BEGIN { FS="\n" }

# define the common suffix as a global variable
suffix = ""

# iterate over each line (i.e. word) in the input file
{
    # set the current word as the value of the first field
    word = $1

    # if the current word is the first word in the file,
    # set the common suffix to be the entire word
    if (NR == 1) {
        suffix = word
    }
    # otherwise, find the longest common suffix between the current word
    # and the current common suffix
    else {
        # initialize the current common suffix to be the entire current word
        curr_suffix = word

        # iterate over each character in the current common suffix,
        # starting from the end of the suffix
        for (i = length(suffix); i > 0; i--) {
            # if the current character of the current common suffix
            # is not present at the corresponding position in the current word,
            # remove the character from the current common suffix
            if (substr(suffix, i, 1) != substr(word, i, 1)) {
                curr_suffix = substr(curr_suffix, 1, i - 1)
            }
        }

        # set the common suffix to be the current common suffix
        suffix = curr_suffix
    }
}

# print the common suffix
END { print suffix }
