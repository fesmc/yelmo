"""Shell completion scripts for the ``yelmo-config`` command.

Emitted by ``yelmo-config completion {bash,zsh}``. They are dependency-free and
use two hidden helper commands (``_complete groups`` / ``_complete params``) to
offer dynamic group / ``group.name`` completion when run inside a checkout.
"""

from __future__ import annotations

BASH = r"""
# bash completion for yelmo-config
# install:  yelmo-config completion bash > /etc/bash_completion.d/yelmo-config
#       or: source <(yelmo-config completion bash)
_yelmo_config_complete() {
    local cur prev sub i
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    local subcommands="groups list show search write diff check update completion"
    local global_opts="--defaults --src --no-color --format --version --help"

    if [[ "$prev" == "--format" ]]; then
        COMPREPLY=( $(compgen -W "text json compact" -- "$cur") ); return
    fi
    if [[ "$prev" == "--defaults" || "$prev" == "--src" || "$prev" == "--root" || "$prev" == "-o" || "$prev" == "--output" ]]; then
        COMPREPLY=( $(compgen -f -- "$cur") ); return
    fi

    sub=""
    for ((i=1; i < COMP_CWORD; i++)); do
        case "${COMP_WORDS[i]}" in
            -*) ;;
            *) sub="${COMP_WORDS[i]}"; break;;
        esac
    done

    if [[ -z "$sub" ]]; then
        COMPREPLY=( $(compgen -W "$subcommands $global_opts" -- "$cur") ); return
    fi

    case "$sub" in
        list)
            COMPREPLY=( $(compgen -W "$(yelmo-config _complete groups 2>/dev/null) --format" -- "$cur") );;
        show)
            COMPREPLY=( $(compgen -W "$(yelmo-config _complete params 2>/dev/null) --format" -- "$cur") );;
        search)
            COMPREPLY=( $(compgen -W "--format" -- "$cur") );;
        write)
            COMPREPLY=( $(compgen -W "-o --output --no-align --format" -- "$cur") $(compgen -f -- "$cur") );;
        diff)
            COMPREPLY=( $(compgen -W "--raw --exit-code --format" -- "$cur") $(compgen -f -- "$cur") );;
        check)
            COMPREPLY=( $(compgen -W "--files --root --format" -- "$cur") $(compgen -f -- "$cur") );;
        update)
            COMPREPLY=( $(compgen -W "--dry-run" -- "$cur") );;
        completion)
            COMPREPLY=( $(compgen -W "bash zsh" -- "$cur") );;
        *)
            COMPREPLY=( $(compgen -W "$global_opts" -- "$cur") );;
    esac
}
complete -F _yelmo_config_complete yelmo-config
"""

ZSH = r"""
#compdef yelmo-config
# zsh completion for yelmo-config
# install:  yelmo-config completion zsh > "${fpath[1]}/_yelmo-config"
#       or: source <(yelmo-config completion zsh)
_yelmo_config() {
    local -a subcommands
    subcommands=(
        'groups:list canonical groups and component mapping'
        'list:list parameters'
        'show:show one parameter (group.name)'
        'search:search names and docs'
        'write:write a complete par file'
        'diff:compare par files'
        'check:validate a par file'
        'update:self-update yelmo-config'
        'completion:emit a shell completion script'
    )
    local state line
    _arguments -C \
        '--defaults[path to yelmo_defaults.nml]:file:_files' \
        '--src[path to src dir]:dir:_files -/' \
        '--no-color[disable coloured output]' \
        '--format[output format]:format:(text json compact)' \
        '--version[show version]' \
        '1:command:->cmd' \
        '*::arg:->args'
    case $state in
        cmd) _describe 'command' subcommands ;;
        args)
            case $line[1] in
                list)   compadd -- ${(f)"$(yelmo-config _complete groups 2>/dev/null)"} ;;
                show)   compadd -- ${(f)"$(yelmo-config _complete params 2>/dev/null)"} ;;
                write)  _arguments '-o[output]:file:_files' '--output[output]:file:_files' '--no-align' '--format:format:(text json compact)' '*:file:_files' ;;
                diff)   _arguments '--raw' '--exit-code' '--format:format:(text json compact)' '*:file:_files' ;;
                check)  _arguments '--files' '--root:dir:_files -/' '--format:format:(text json compact)' '*:file:_files' ;;
                update) _arguments '--dry-run' '1:ref:' ;;
                completion) compadd bash zsh ;;
            esac ;;
    esac
}
_yelmo_config "$@"
"""


def script(shell: str) -> str:
    s = {"bash": BASH, "zsh": ZSH}.get(shell)
    if s is None:
        raise ValueError(f"unsupported shell '{shell}' (choose bash or zsh)")
    return s.lstrip("\n")
