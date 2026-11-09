# Discussion 9: Debugging

Work with your group to investigate two programs. Record the evidence that leads you to a change and check the result after making it.

## Setup

Run all shell commands below from this directory (the repository root, containing instructions.md and Makefile).

```bash
make environment
```

The Python environment includes the visualization dependencies. The C++ debugging activity additionally requires `make`, a C++17 compiler (`g++`), and GDB on a Linux machine or Linux environment. These system tools are not installed by `make environment`.

If you are not using Linux or WSL, you can run these exercises on GitHub Codespaces.
If you are running on GitHub Codespaces, you will have to run **conda init**, then open a new terminal before you can activate the environment.

Activate the environment with:

```
conda activate discussion-9
```

Check before starting:

```bash
python --version
g++ --version
gdb --version
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

Remove all breakpoints for the next part.

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
3. For each error, record the breakpoint location, the commands that provided useful evidence, and what that evidence showed.
4. Make a change and rerun. You can remove the breakpoint to rerun if you would like. Continue until the program produces an image.

For example, the first error occurs at `data = text.split("\n")` (line 28 in the supplied file). After you observe this error, insert `breakpoint()` immediately above that statement, inside the `with` block. Run the program again and examine the variables and data types in the debugger. Try `p type(text)` and consider how to convert the value to the appropriate data type. Line numbers will change as you edit the file.

How would your approach change if the failing operation were inside a library you did not want to edit? Identify a breakpoint in your own code that would let you inspect the values passed to that library.

### Check the image

The output is `caffeine.png` in your current directory. Compare it with `python_debugging/caffeine_reference.png`. For the local fixture, the graph should contain **14 non-hydrogen atoms and 15 bonds**. Check element labels, connectivity, colors and bond widths; an exact pixel match is not required.

The script calculates a layout. Is that layout actually used? Make any changes needed to ensure that the graph layout is used. What kind of tool would have caught this mistake (hint: we learned about it in Week 1)?

Install the tool you identified and run a check with it on the Python file. Does it identify any other issues?

## 2. GDB: find the closest pair of atoms

`gdb_debugging/closest_pair.cpp` reads an XYZ file and reports the closest pair and their distance in angstrom. Atom indices in the output are zero-based. The program compiles but terminates during the calculation.

First read the code. Identify the file-reading function, the distance calculation, and the function that searches over pairs.

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

Enter the following commands at the `(gdb)` prompt, rather than in your shell. A command reference table appears below, just before the breakpoint investigation.

First, start the program:

```text
run
```

If you see the following prompt, answer **n** and press Enter. Your executable already includes debugging information from `-g`; downloading system-library debugging information is not required for this exercise.

```text
This GDB supports auto-downloading debuginfo from the following URLs:
  <https://debuginfod.ubuntu.com>
Enable debuginfod for this session? (y or [n])
```

After `run`, the program stops with an error. You may see `SIGABRT` and a warning about a missing system source file such as `pthread_kill.c`. Continue with your investigation; you do not need that file.

Next, type `bt` (short for `backtrace`) to show the stack of function calls that led to where the program stopped:

```text
bt
```

Look for `find_closest_pair` in this output. Each entry is a **stack frame**, representing a function call. Its frame number is the number immediately after `#` at the start of the entry. For example, if you see:

```text
#10 find_closest_pair (...) at gdb_debugging/closest_pair.cpp:43
```

select that frame with `frame 10` (without the `#`). Your frame number may differ; use the number shown in your own backtrace. Then inspect its source and variables:

```text
list
info locals
print i
print j
print n_atoms
```

Next, use the exception workflow from the lecture to stop when the exception is thrown, before the program aborts. Stay in the same GDB session and enter these commands at the `(gdb)` prompt:

```text
catch throw
run
bt
```

`catch throw` tells GDB to pause whenever a C++ exception is thrown, before the program handles it or terminates because it is unhandled. This is an exception catchpoint: you do not need to know the failing source line in advance. Here, it lets you inspect the call stack and variables when the bounds-checking exception is thrown, rather than after the program aborts.

`run` restarts your program with the catchpoint active. If GDB asks whether to start it from the beginning, answer **y**. When the catchpoint stops execution, `bt` shows the call stack at the throw.

Again find `find_closest_pair` in the new backtrace and select its frame using the number after `#`. The frame number can differ from the earlier backtrace. Inspect `i`, `j`, and `n_atoms` again. You do not need to repair C++ library internals.

Use this command table as a reference as you investigate with breakpoints:

| Command | Action |
| --- | --- |
| `run` | Start or restart the executable. |
| `catch throw` | Pause when a C++ exception is thrown. |
| `bt` | Show the backtrace (the stack of function calls). |
| `frame NUMBER` | Select the stack frame labeled `#NUMBER` in the backtrace. |
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

## Report out and submission

Submit your corrected Python visualization and C++ closest-pair code, and a short group README with:

- One useful piece of debugger evidence for each bug, the resulting fix, and how you checked it.
- Your visualization checks and closest-pair results.
- The tool you used to check the Python file, what it reported, and the changes you made.

Remove Python debugging breakpoints before submitting. On your Report Out slide, summarize one debugging discovery and the evidence that explains why your change worked.

## Documentation

- [Python debugger](https://docs.python.org/3.15/library/pdb.html)
- [GDB backtraces](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Backtrace.html)
- [GDB exception catchpoints](https://sourceware.org/gdb/current/onlinedocs/gdb.html/Set-Catchpoints.html)
- [NetworkX drawing](https://networkx.org/documentation/stable/reference/generated/networkx.drawing.nx_pylab.draw_networkx.html)
- [urllib.request](https://docs.python.org/3.15/library/urllib.request.html)
