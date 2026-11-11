_varigraph_completion()
{
    COMPREPLY=()

    if [[ $COMP_CWORD -eq 1 ]]; then
        COMPREPLY=($(compgen -W 'construct genotype' -- "${COMP_WORDS[COMP_CWORD]}"))
    fi
}

complete -o default -F _varigraph_completion varigraph varigraph-gpu
