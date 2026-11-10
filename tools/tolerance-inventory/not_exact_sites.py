#!/usr/bin/env python3
"""List non-test source sites that name or silently consume `Outcome::NotExact`.

`explicit` rows name `NotExact`. `implicit` rows match an `Outcome` without a
`NotExact` arm, so a wildcard, binding or `else` arm receives it. The scan is a
textual heuristic: it misses method-style consumers such as `.is_coincident()`
and reports a few `matches!` tests that are not on an `Outcome`. Use it to find
sites that docs/NOT-EXACT-INVENTORY.md does not yet cover, not as a proof.
"""
import argparse
from pathlib import Path
import re

SKIPPED = ('/tests/', '/benches/', '/examples/')
FUNCTION = re.compile(r'\s*(pub(\([^)]*\))?\s+)?(const\s+)?(async\s+)?fn\s+(\w+)')
OUTCOME_PATTERN = re.compile(r'\b(Outcome::|Proven\(|ProvenCoincident\b)')
CATCH_ALL = re.compile(r'(_|[a-z_][a-z0-9_]*)(\s+if\s.*)?', re.S)
LET_FORM = re.compile(r'\b(if|while)\s+let\s+(catallax_numbers::)?Outcome::'
                      r'|\blet\s+(catallax_numbers::)?Outcome::')


def sources(root):
    for path in sorted((root / 'crates').rglob('*.rs')):
        rel = path.relative_to(root).as_posix()
        if any(part in rel for part in SKIPPED) or rel.endswith('tests.rs'):
            continue
        lines = path.read_text().split('\n')
        # Inline test modules conventionally close the file.
        for index, line in enumerate(lines):
            if (line.strip() == '#[cfg(test)]' and index + 1 < len(lines)
                    and re.match(r'\s*mod \w+\s*\{', lines[index + 1])):
                lines = lines[:index]
                break
        yield rel, '\n'.join(lines)


def blank_comments_and_strings(text):
    out = list(text)
    i, n = 0, len(text)
    while i < n:
        if text.startswith('//', i):
            end = text.find('\n', i)
            end = n if end < 0 else end
            out[i:end] = ' ' * (end - i)
            i = end
        elif text[i] == '"':
            j = i + 1
            while j < n and text[j] != '"':
                j += 2 if text[j] == '\\' else 1
            for k in range(i + 1, min(j, n)):
                if text[k] != '\n':
                    out[k] = ' '
            i = j + 1
        else:
            i += 1
    return ''.join(out)


def line_of(text, offset):
    return text.count('\n', 0, offset) + 1


def match_arms(text, start):
    """Top-level `pattern => body` pairs of the match block opening at `start`."""
    depth, end = 0, start
    while end < len(text):
        depth += {'{': 1, '}': -1}.get(text[end], 0)
        if depth == 0:
            break
        end += 1
    body, arms, depth, k, pattern_start = text[start + 1:end], [], 0, 0, 0
    while k < len(body):
        c = body[k]
        if c in '([{':
            depth += 1
        elif c in ')]}':
            depth -= 1
        elif depth == 0 and body.startswith('=>', k):
            q, inner = k + 2, 0
            while q < len(body):
                cc = body[q]
                if cc in '([{':
                    inner += 1
                elif cc in ')]}':
                    inner -= 1
                    if inner == 0 and cc == '}' and body[k + 2:q + 1].strip().startswith('{'):
                        q += 1
                        break
                elif cc == ',' and inner == 0:
                    break
                q += 1
            arms.append((body[pattern_start:k].strip(), body[k + 2:q].strip()))
            k = pattern_start = q + 1
            continue
        k += 1
    return arms


def explicit_sites(rel, text):
    function = '?'
    for number, line in enumerate(text.split('\n'), 1):
        found = FUNCTION.match(line)
        if found:
            function = found.group(5)
        if 'NotExact' in line:
            yield rel, number, 'explicit', function, line.strip()


def implicit_sites(rel, raw):
    text = blank_comments_and_strings(raw)
    for found in re.finditer(r'\bmatch\b', text):
        i, depth = found.end(), 0
        while i < len(text):
            c = text[i]
            if c in '([':
                depth += 1
            elif c in ')]':
                depth -= 1
            elif depth == 0 and c in '{;':
                break
            i += 1
        if i >= len(text) or text[i] != '{':
            continue
        arms = match_arms(text, i)
        patterns = [pattern for pattern, _ in arms]
        if (not any(OUTCOME_PATTERN.search(p) for p in patterns)
                or any('NotExact' in p for p in patterns)):
            continue
        catch = [arm for arm in arms if CATCH_ALL.fullmatch(arm[0])]
        if catch:
            body = re.sub(r'\s+', ' ', catch[0][1])[:100]
            yield rel, line_of(text, found.start()), 'implicit', 'match', body
    for found in re.finditer(r'matches!\s*\(', text):
        segment = text[found.end():found.end() + 200]
        if OUTCOME_PATTERN.search(segment) and 'NotExact' not in segment:
            yield (rel, line_of(text, found.start()), 'implicit', 'matches!',
                   re.sub(r'\s+', ' ', segment)[:100])
    for found in LET_FORM.finditer(text):
        yield (rel, line_of(text, found.start()), 'implicit', 'let',
               re.sub(r'\s+', ' ', text[found.start():found.start() + 100]))


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    rows = []
    for rel, text in sources(args.root):
        rows.extend(explicit_sites(rel, text))
        rows.extend(implicit_sites(rel, text))
    for rel, number, kind, where, excerpt in sorted(rows):
        print(f'{rel}:{number}\t{kind}\t{where}\t{excerpt}')


if __name__ == '__main__':
    main()
