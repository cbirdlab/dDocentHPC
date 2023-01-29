BEGIN {
    FS = "\n"
}

{
    words[NR] = $0
}

END {
    common_suffix = ""
    min_length = length(words[1])

    for (i = 2; i <= NR; i++) {
        if (length(words[i]) < min_length) {
            min_length = length(words[i])
        }
    }

    for (i = min_length; i >= 1; i--) {
        current_char = substr(words[1], i, 1)
        for (j = 2; j <= NR; j++) {
            if (substr(words[j], i, 1) != current_char) {
                print common_suffix
                exit
            }
        }
        common_suffix = current_char common_suffix
    }

    print common_suffix
}
