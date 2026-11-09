# Discussion 9: Debugging and Profiling

Work with your group to investigate three programs. Record the evidence that leads you to a change and check the result after making it.

## Setup

Run all shell commands below from this directory (the repository root, containing instructions.md and Makefile).

```bash
make environment
conda activate discussion-9
```

The Python environment includes the visualization dependencies. The C++ debugging and profiling activities additionally require `make`, a C++17 compiler (`g++`), GDB, and GNU gprof on a Linux machine or Linux environment. These system tools are not installed by `make environment`.

If you are not using Linux or WSL, you can run these exercises on GitHub Codespaces.

Check before starting:

```bash
python --version
g++ --version
gdb --version
gprof --version
make --version
```

## 1. Python debugging: molecular visualization

This script retrieves a molecule by name from PubChem and draws its molecular graph using NetworkX. This is similar to your molecular graph work in Problem Set 3. The supplied script contains errors.

### Read before running

Read `python_debugging/main.py` and trace what happens when it is invoked with `caffeine` as its first argument. Summarize:

1. What does each function receive and return? In what order are the functions called?
2. Where is a context manager used? What resource does it manage? Are any decorators present?
3. How does the script obtain command-line arguments? What does its optional second argument do?
4. What does a graph node represent, and what identifies an individual atom?

### Explore the debugger

| Command | Action |
| --- | --- |
| `n` | Execute the next line, stepping over function calls. |
| `s` | Step into a Python function call. |
| `c` | Continue to the next breakpoint or the end. |
| `p expression` | Evaluate and print an expression. |
| `l` | List nearby source code. |
| `w` | Show the call stack. |
| `u` / `d` | Move to the caller / a newer stack frame. |
| `q` | Quit the debugger. |

Place `breakpoint()` near the top of the module, before the imports. Run:

```bash
python python_debugging/main.py
```

Use `l` to identify the line that will execute next. Step with `n`, then quit. Move the breakpoint to just inside `if __name__ == "__main__":` and repeat. What code has already executed when you reach this breakpoint? What happens when execution reaches the argument access without a molecule name?

### Find and repair the failures

Now supply a molecule name:

```bash
python python_debugging/main.py caffeine
```

If PubChem is unavailable, use the supplied local structure. This exercises the same parsing and graph-building code, including the intended bugs:

```bash
python python_debugging/main.py caffeine data/caffeine.sdf
```

For each failure:

1. Read the traceback and form a hypothesis.
2. Place a breakpoint before the relevant operation. Use `p`, `type()`, `dir()`, and stepping to examine the inputs and intermediate values.
3. Record the breakpoint location, the commands that provided useful evidence, and what that evidence showed.
4. Make a change and rerun. Continue until the program produces an image.

How would your approach change if the failing operation were inside a library you did not want to edit? Identify a breakpoint in your own code that would let you inspect the values passed to that library.

### Check the image

The output is `caffeine.png` in your current directory. Compare it with `python_debugging/caffeine_reference.png`. For the local fixture, the graph should contain **14 non-hydrogen atoms and 15 bonds**. Check element labels, connectivity, colors and bond widths; an exact pixel match is not required.

The script calculates a layout. Is that layout actually used? Also trace how the bond-width list is matched to the graph's edge order. Make any changes needed to ensure that each bond gets the width associated with its own order.

Explain why producing an image without raising an exception is insufficient evidence of correctness.

### Code improvement

Add handling for a molecule that PubChem cannot find. Make sure a network failure is not automatically reported as “molecule not found.” Decide how the program should respond to a missing command-line argument. Improve documentation where useful. Remove debugging breakpoints before submitting.

Optional: support an element absent from the current color map with a sensible default.

## 2. GDB: find the closest pair of atoms

`gdb_debugging/closest_pair.cpp` reads an XYZ file and reports the closest pair and their distance in angstrom. Atom indices in the output are zero-based. The program compiles but terminates during the calculation.

First read the code. Identify the file-reading function, the distance calculation, and the function that searches over pairs. Why should an atom not be compared with itself?

### Build and reproduce

```bash
make gdb
./gdb_debugging/closest_pair data/check.xyz
```

The Makefile uses `-g` to include debugging information and `-O0` to disable compiler optimization for this exercise. Record the termination message.

### Investigate in GDB

```bash
gdb --args ./gdb_debugging/closest_pair data/check.xyz
```

At the GDB prompt:

```text
run
bt
```

Look for `find_closest_pair` in the backtrace. Its frame number will vary. Select that frame using `frame NUMBER`, substituting its actual number, then inspect:

```text
list
info locals
print i
print j
print n_atoms
```

Next, restart GDB and stop at the exception before termination, using the exception workflow from the lecture:

```text
catch throw
run
bt
```

Again select the frame for `find_closest_pair`. You do not need to repair C++ library internals.

| Command | Action |
| --- | --- |
| `run` | Start or restart the executable. |
| `bt` | Show the backtrace. |
| `frame NUMBER` | Select a stack frame. |
| `up` / `down` | Move between caller and callee frames. |
| `list` | Show nearby source. |
| `print expression` | Inspect an expression in the selected frame. |
| `info locals` | Show local variables. |
| `break find_closest_pair` | Stop on entry to this function. |
| `break closest_pair.cpp:LINE` | Set a breakpoint at a source line. |
| `next` / `step` | Step over / into calls. |
| `continue` | Resume execution. |
| `info breakpoints` | List breakpoints and catchpoints. |
| `delete NUMBER` | Remove a breakpoint or catchpoint. |
| `quit` | Exit GDB. |

Choose a breakpoint before the failing access, restart, and step through the loop to see how the invalid state is reached. Optional: use `break closest_pair.cpp:LINE if CONDITION` to stop only on the iteration you want to examine.

Answer:

1. Which operation throws the exception? Which line of your code causes it?
2. What are the relevant index values and the number of atoms at that point?
3. What are the valid indices? How did the loop reach an invalid one?
4. Why did you need to select a different stack frame to inspect these variables?

### Repair and verify

Edit the source, exit GDB, and rebuild with `make gdb`. Start a new GDB session or run the rebuilt executable normally. Editing a `.cpp` file does not update an already-built executable.

For `data/check.xyz`, the correct closest pair is **1 and 3**, at **2.000000 angstrom**. Verify this by hand from the coordinates as well as by running the program.

Then run:

```bash
./gdb_debugging/closest_pair data/caffeine.xyz
```

The caffeine XYZ contains hydrogens; the visualization omits them by default. Explain why the closest pair need not correspond to the shortest-looking edge in the drawn graph.

## 3. C++ profiling with gprof: minimum distances across structures

`profiling/profile_distances.cpp` calculates the minimum interatomic distance for each generated configuration. The program gives correct results but does unnecessary work. Inputs are synthetic points on a perturbed lattice in angstrom, not an MD trajectory. It generates them once using a fixed seed and repeats the analysis to obtain a measurable runtime.

### Predict, then profile

Read the program and predict its most expensive function. Explain your prediction before profiling. Follow the workflow from lecture: compile with `-pg`, run, then inspect `gmon.out` using `gprof`.

From the directory containing the Makefile:

```bash
make profile
./profiling/distances_before 120 12 10 profiling/before.txt
gprof ./profiling/distances_before gmon.out > profiling/before_profile.txt
mv gmon.out profiling/before_gmon.out
```

The three numeric arguments are **atoms per structure, number of structures, and number of analysis repetitions**. This build uses `-std=c++17 -g -O0 -pg`. The one-step compile command applies `-pg` during both compilation and linking. Run the program to normal completion so it writes `gmon.out`.

Read `profiling/before_profile.txt`. Standard-library routines may appear prominently in this unoptimized build. Use the call graph to connect their work to the application functions that call them.

| Output | What to inspect |
| --- | --- |
| Flat profile: `% time`, `self seconds` | Where sampled execution time is spent in the function itself. |
| Flat profile: `calls` | How often each instrumented function was called. |
| Flat profile: `cumulative seconds` | A running total down the table, **not** the function's time including its callees. |
| Call graph | Callers, descendants, and the separate self and children times. |

Answer:

1. Does the profile support your prediction? Identify an expensive path through your own code.
2. How often is `pairwise_distances` called? Explain the count using the three input arguments.
3. Which function calls it? Do its inputs change between all of these calls?
4. Why can `minimum_distance` have little self time but substantial time in its children?
5. How does the meaning of gprof's flat-profile “cumulative seconds” differ from Python cProfile's `cumtime`?

### Make one targeted improvement

Copy `profiling/profile_distances.cpp` to `profiling/profile_distances_optimized.cpp`. Keep the original as your baseline. Use your measurements to eliminate unnecessary computation while preserving the results. Do not change the input generator, distance formula, or command-line arguments.

```bash
g++ -std=c++17 -g -O0 -pg profiling/profile_distances_optimized.cpp -o profiling/distances_after
./profiling/distances_after 120 12 10 profiling/after.txt
gprof ./profiling/distances_after gmon.out > profiling/after_profile.txt
mv gmon.out profiling/after_gmon.out
diff profiling/before.txt profiling/after.txt
```

No output from `diff` indicates that all reported distances agree. Also compare the printed checksums. Compare the distance-matrix call counts in the two profiles, and explain the change.

Each run writes `gmon.out` in the current directory. Keep the before and after files separate, and always analyze a profile with the exact executable that produced it. If a short run reports little or no sampled time, that does not mean the functions take no time. Increase the repetition count for **both** versions and repeat the comparison. Call counts can still be informative when timing samples are sparse.

### Measure performance without instrumentation

Profiling adds overhead. To judge performance, compile both versions with the same optimization settings and without `-pg`, as discussed in lecture:

```bash
g++ -std=c++17 -O2 profiling/profile_distances.cpp -o profiling/before_release
g++ -std=c++17 -O2 profiling/profile_distances_optimized.cpp -o profiling/after_release
```

Run each version three times with the same inputs:

```bash
./profiling/before_release 120 12 100 profiling/before_release.txt
./profiling/after_release 120 12 100 profiling/after_release.txt
```

Repeat that pair of commands twice more. Record the median “Analysis seconds” for each version and calculate speedup as `before / after`. The timer covers repeated analysis and checksum accumulation; it excludes coordinate generation and file output. Check the output files with `diff` and compare checksums again.

Report the timings, input arguments, compiler flags and call counts. Do not compare a baseline `-O0 -pg` time against an improved `-O2` time to claim an algorithmic speedup.

Optional: predict how doubling the atom count changes the amount of work in each version, then investigate using 240 atoms. Keep all other settings the same in the comparison.

## Report out and submission

Submit your corrected visualization and C++ code, your optimized C++ profiling source, and a short group README with:

- One useful piece of debugger evidence for each bug, the resulting fix, and how you checked it.
- The profiler's relevant call counts before and after the change.
- Matching-output evidence, benchmark times and measured speedup.
- Your error-handling decisions for the visualization script.

On your Report Out slide, summarize one debugging discovery and the profiling result. Include enough evidence to explain why your change worked; you do not need a transcript of every command.

## Documentation

- [Python debugger](https://docs.python.org/3.11/library/pdb.html)
- [GNU gprof manual](https://sourceware.org/binutils/docs/gprof/)
- [GDB backtraces](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Backtrace.html)
- [GDB exception catchpoints](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Set-Catchpoints.html)
- [NetworkX drawing](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.nx_pylab.draw_networkx.html)
- [urllib.request](https://docs.python.org/3.11/library/urllib.request.html)
