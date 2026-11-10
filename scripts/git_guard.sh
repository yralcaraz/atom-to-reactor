#!/usr/bin/env bash
# Authorship and raw-data guard of the repository, called by the hooks in scripts/hooks/
# (enable once per clone: git config core.hooksPath scripts/hooks).
#
#   scripts/git_guard.sh pre-commit              staged raw files; author and committer identity
#   scripts/git_guard.sh commit-msg <file>       Co-Authored-By / claude / anthropic lines in the message
#   scripts/git_guard.sh pre-push <remote>       every commit to be pushed (refs on stdin, as git gives them)
#   scripts/git_guard.sh range <rev-range>       the same checks on a range, by hand (e.g. origin/M1-steps..HEAD)
#
# Raw data (never in git): the NMR machine files. Source: Y. Alcaraz Galván

set -u
NAME='yralcaraz'
EMAIL='yerayalc@gmail.com'
RAW='(\.(jdf|jdx|jcamp|fid|ser|mnova)$|(^|/)(fid|ser|acqu|acqus|procpar)$)'
MSG='co-authored-by|claude|anthropic'
FAIL=0

refuse() { echo "git_guard: REFUSED: $*" >&2; FAIL=1; }

check_ident() {   # $1 role, $2 "name <email>"
    [[ "$2" == "$NAME <$EMAIL>" ]] || refuse "$1 is '$2', not '$NAME <$EMAIL>'"
}

check_range() {   # every commit of a rev-list argument set
    local c
    for c in $(git rev-list "$@"); do
        check_ident "author of ${c:0:7}" "$(git show -s --format='%an <%ae>' "$c")"
        check_ident "committer of ${c:0:7}" "$(git show -s --format='%cn <%ce>' "$c")"
        git show -s --format=%B "$c" | grep -iqE "$MSG" && refuse "message of ${c:0:7} has a line with claude/anthropic/co-authored"
        raw=$(git diff-tree --no-commit-id --root -r --name-only --diff-filter=ACMR "$c" | grep -iE "$RAW")
        [[ -z "$raw" ]] || refuse "${c:0:7} adds raw files: $raw"
    done
}

case "${1:-}" in
    pre-commit)
        check_ident author "$(git var GIT_AUTHOR_IDENT | sed 's/ [0-9]* [-+][0-9]*$//')"
        check_ident committer "$(git var GIT_COMMITTER_IDENT | sed 's/ [0-9]* [-+][0-9]*$//')"
        raw=$(git diff --cached --name-only --diff-filter=ACMR | grep -iE "$RAW")
        [[ -z "$raw" ]] || refuse "raw files staged: $raw"
        ;;
    commit-msg)
        grep -v '^#' "$2" | grep -iqE "$MSG" && refuse "the message has a line with claude/anthropic/co-authored"
        ;;
    pre-push)
        zero=0000000000000000000000000000000000000000
        while read -r lref lsha rref rsha; do
            [[ "$lsha" == "$zero" ]] && continue                      # deletion
            [[ "$lref" == refs/tags/backup/* || "$rref" == refs/tags/backup/* ]] && { refuse "tag $lref is a backup tag"; continue; }
            if [[ "$rsha" == "$zero" ]]; then check_range "$lsha" --not --remotes   # new branch: what no remote has yet
            else check_range "$rsha..$lsha"; fi
        done
        ;;
    range)
        shift; check_range "$@"
        ;;
    *)
        echo "usage: $0 pre-commit | commit-msg <file> | pre-push <remote> | range <rev-range>" >&2; exit 2
        ;;
esac
[[ $FAIL == 0 ]] || exit 1
