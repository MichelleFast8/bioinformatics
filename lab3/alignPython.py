# Michelle Fast, Andrew Holmes
# Bioinformatics Lab 3

import itertools
import numpy as np
import sys


def calculate_score(src_word, matching_word):
    score = 0
    match, mismatch = 5, -3
    for (a, b) in zip(src_word, matching_word):
        if (a==b):
            score += match
        else:
            score += mismatch
    return score


def main():

    # Parse -g flag
    gaps = False
    if sys.argv[-1] == '-g':
        print("Allowing gaps")
        gaps = True

    # Ask user for Source sequence
    src = input("SOURCE: ")

    # Create words
    w = 3
    src_words = [src[i] + src[i+1] + src[i+2] for i in range(len(src)-w+1)]
    matching_words = list(itertools.product("ATCG", repeat=w))
    matching_words = ["".join(word) for word in matching_words]

    # Compare each src word to all matching words, store in matrix
    scored_words = np.zeros((len(matching_words), len(src_words))) # Rows = matching words, Cols = src words
    for j, src_word in enumerate(src_words):
        for i, matching_word in enumerate(matching_words):
            score = calculate_score(src_word, matching_word)
            scored_words[i][j] = score

    # Take only the matching_words that align best with src_words
    p = 0.3
    indices = (scored_words > 0).sum(axis=1).argsort()[::-1][:int(len(matching_words) * p)]

    # Store top matching words along with corresponding src words with positive score
    top_matches = []
    for index in indices:
        m_word = matching_words[index]
        s_words = np.take(src_words, [idx for idx, score in enumerate(scored_words[index]) if score > 0])
        top_matches.append((m_word, s_words))

    # Compare to database!
    VERY_BEST = {
            "seq_id": "",
            "db_sequence": "",
            "src": "",
            "score": -np.inf,
            }

    db = open("seq.txt")
    for line in db:

        # Clean up database entry
        entry_id = line.split(" ")[0]
        seq = line.split(" ")[1][:-1]

        # Loop through matching words
        for match in top_matches:

            # See if matching word is in db entry
            db_sequence_index = seq.find(match[0])
            if (db_sequence_index == -1):
                continue

            # Compare each source word with db entry
            for i, src_word in enumerate(match[1]):

                # Index of src_word in source
                query_index = src.find(src_word)

                # If no gaps
                if not gaps:
                    padded_seq = len(src) * ' ' + seq + len(src) * ' '
                    db_sequence_index = db_sequence_index + len(src)
                    padded_src = (db_sequence_index - query_index) * ' ' + src
                    padded_seq_slice = padded_seq[db_sequence_index - query_index: db_sequence_index - query_index + len(src)]

                    score = calculate_score(src, padded_seq_slice)

                    # Update overall best match
                    if score >= VERY_BEST['score']:
                        VERY_BEST['score'] = score
                        VERY_BEST['db_sequence'] = padded_seq
                        VERY_BEST['src'] = padded_src
                        VERY_BEST['seq_id'] = entry_id

                else:

                    # Score alignment
                    score = calculate_score(match[1][i], match[0])

                    # Extend alignment right
                    query_right, db_sequence_right = query_index + w, db_sequence_index + w
                    while(query_right < len(src) and db_sequence_right < len(seq)):

                        # Compare extension
                        if src[query_right] == seq[db_sequence_right]:
                            score += 5
                            query_right += 1
                            db_sequence_right += 1
                            src_word += src[query_right]

                        else:
                            # gaps here?
                            break

                    # Extend alignment left
                    query_left, db_sequence_left = query_index - 1, db_sequence_index - 1
                    while(query_left > 0 and db_sequence_left > 0):

                        # Compare extension
                        if src_word[query_left] == seq[db_sequence_left]:
                            score += 5
                            query_left -= 1
                            db_sequence_left -= 1
                            src_word = src[query_left] + src_word

                        else:
                            # gaps here?
                            break


                    # Update overall best match
                    if score >= VERY_BEST['score']:
                        VERY_BEST['score'] = score
                        VERY_BEST['db_sequence'] = seq
                        VERY_BEST['src_index'] = db_sequence_left
    breakpoint()
    # print(VERY_BEST)



if __name__ == "__main__":
    main()
