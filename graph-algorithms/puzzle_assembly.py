# python3
"""Problem 1: Puzzle Assembly — backtracking, row by row."""
import sys

BLACK = "black"


def parse(line):
    return tuple(line.strip().replace(" ", "").strip("()").split(","))  # (up, left, down, right)


def solve(pieces):
    n = int(round(len(pieces) ** 0.5))
    grid = [[None] * n for _ in range(n)]
    used = [False] * len(pieces)

    def fits(p, i, j):
        up, left, down, right = p
        # Border edges must be black, inner edges must not be.
        if (up == BLACK) != (i == 0): return False
        if (left == BLACK) != (j == 0): return False
        if (down == BLACK) != (i == n - 1): return False
        if (right == BLACK) != (j == n - 1): return False
        # Match neighbours already placed.
        if i > 0 and grid[i - 1][j][2] != up: return False
        if j > 0 and grid[i][j - 1][3] != left: return False
        return True

    def place(cell):
        if cell == n * n:
            return True
        i, j = divmod(cell, n)
        for idx, p in enumerate(pieces):
            if not used[idx] and fits(p, i, j):
                used[idx] = True
                grid[i][j] = p
                if place(cell + 1):
                    return True
                used[idx] = False
        grid[i][j] = None
        return False

    place(0)
    return grid


def main():
    pieces = [parse(l) for l in sys.stdin if l.strip()]
    grid = solve(pieces)
    out = []
    for row in grid:
        out.append(";".join("(" + ",".join(p) + ")" for p in row))
    print("\n".join(out))


if __name__ == "__main__":
    main()
