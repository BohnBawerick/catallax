set -u
echo "== ls-remote"; git ls-remote https://github.com/ZIGRAWBOHNBAWERICK/cad-test-corpus.git HEAD; echo "ls-remote exit $?"
d=$(mktemp -d)
echo "== clone (no checkout, blobless)"; git clone --depth 1 --filter=blob:none --no-checkout https://github.com/ZIGRAWBOHNBAWERICK/cad-test-corpus.git "$d/corpus"; rc=$?; echo "clone exit $rc"
[ $rc = 0 ] && git -C "$d/corpus" log -1 --format='%H %ci'
exit $rc
