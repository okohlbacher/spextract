#!/usr/bin/env bash
# Adversarial review across THREE independent CLIs: codex, kimi, vibe.
#
# Standing instruction: every adversarial review in this project runs all three. Independent
# models catch different things, and this project's failures have consistently been ones the
# author could not see in their own work -- eleven data-handling errors in one session, each
# obvious in hindsight.
#
# Invocation notes, learned the hard way. A first attempt produced NO findings from either CLI
# and I nearly reported "no issues found":
#   codex  needs -o/--output-last-message FILE. Without it stdout carries only the event
#          transcript (skill loading, tool calls) and the actual answer is never printed.
#          Also needs --skip-git-repo-check when run outside a repo.
#   kimi   -p prints to stdout but truncates on long reasoning; --output-format text is
#          explicit, and the response must be captured whole rather than tailed.
#   vibe   -p is programmatic mode; --output text; needs --workdir to read files.
# All three get a generous timeout: a hostile review of a real codebase is not a quick call.
#
# Usage:  adv_review.sh PROMPT_FILE OUTPUT_DIR [FILES_TO_READ...]
set -u

PROMPT_FILE="${1:?usage: adv_review.sh PROMPT_FILE OUTPUT_DIR [paths...]}"
OUT="${2:?need output dir}"
shift 2
mkdir -p "$OUT"
TIMEOUT="${ADV_TIMEOUT:-900}"

[ -s "$PROMPT_FILE" ] || { echo "ABORT: prompt file empty: $PROMPT_FILE" >&2; exit 2; }

# `timeout` is GNU coreutils and is ABSENT on stock macOS -- the first run of this script died
# with 127 on all three reviewers. The yield check below caught it (that is what it is for),
# but a review harness that silently cannot run is exactly the failure mode being guarded
# against, so resolve the binary explicitly and degrade to no timeout rather than to nothing.
if command -v timeout  >/dev/null 2>&1; then TO="timeout $TIMEOUT"
elif command -v gtimeout >/dev/null 2>&1; then TO="gtimeout $TIMEOUT"
else TO=""; echo "[adv] no timeout(1) available -- running without a time limit" >&2
fi
PROMPT="$(cat "$PROMPT_FILE")"
if [ $# -gt 0 ]; then
  PROMPT="$PROMPT

Read these files before answering:
$(printf '  %s\n' "$@")"
fi

echo "[adv] prompt $(wc -c < "$PROMPT_FILE") bytes -> $OUT"

# --- codex ------------------------------------------------------------------------------------
# -o is mandatory: without it the answer never reaches stdout. `< /dev/null` is mandatory too --
# these run backgrounded, and codex exec reads stdin IN ADDITION to the prompt arg, so an
# inherited never-EOF stdin makes it hang on "Reading additional input from stdin..." forever.
( $TO codex exec --skip-git-repo-check -s read-only \
      -c model_reasoning_effort=high \
      -o "$OUT/codex.md" "$PROMPT" > "$OUT/codex.transcript" 2>&1 < /dev/null
  echo "exit=$?" >> "$OUT/codex.transcript" ) &
P_CODEX=$!

# --- kimi -------------------------------------------------------------------------------------
( $TO kimi -p "$PROMPT" --output-format text > "$OUT/kimi.md" 2>"$OUT/kimi.err" < /dev/null
  echo "exit=$?" >> "$OUT/kimi.err" ) &
P_KIMI=$!

# --- vibe -------------------------------------------------------------------------------------
( $TO vibe -p "$PROMPT" --output text --auto-approve \
      --workdir "$(pwd)" > "$OUT/vibe.md" 2>"$OUT/vibe.err" < /dev/null
  echo "exit=$?" >> "$OUT/vibe.err" ) &
P_VIBE=$!

wait $P_CODEX $P_KIMI $P_VIBE

# --- report which reviewers actually produced findings -----------------------------------------
# A silent CLI must be visible as SILENT, not folded into "no issues found". The first attempt
# at this returned two empty reviews and the absence nearly went unremarked.
echo
echo "=== reviewer yield ==="
STATUS=0
for R in codex kimi vibe; do
  F="$OUT/$R.md"
  N=$( [ -f "$F" ] && wc -l < "$F" || echo 0 )
  B=$( [ -f "$F" ] && wc -c < "$F" || echo 0 )
  if [ "$B" -lt 200 ]; then
    echo "  $R: SILENT (${B} bytes) -- treat as NO REVIEW, not as a clean bill of health"
    [ -s "$OUT/$R.err" ] && sed -n '1,3p' "$OUT/$R.err" | sed 's/^/      /'
    STATUS=1
  else
    echo "  $R: $N lines, $B bytes"
  fi
done
echo
echo "outputs: $OUT/{codex,kimi,vibe}.md"
exit $STATUS
