# Python Debugging  

To begin, create an environment for this discussion by using the Makefile:

```
make environment
```

Then, activate the environment:

```
conda activate discussion-9
```

This discussion contains a Python script which can visualize a molecule retrieved by name from [PubChem](https://pubchem.ncbi.nlm.nih.gov/). 
This task is very similar to what you implemented with molecular graphs in Problem Set 3. 
Below, you will apply the script to caffeine. A reference image for caffeine is in the file `caffeine_reference.png`. 
The script provided for this discussion contains some errors. 

Your task will be to *[use the Python debugger](https://docs.python.org/3/library/pdb.html)* to find and fix the problems. 

To use the script to retrieve caffeine, type the following into the shell:

```
python main.py caffeine
```

You will notice that if you try the command initially, the script will not work. Part of your task is to debug and fix this problem.

Before you debug the code, look over the script and answer the following questions (put summaries of your answers on your Report Out slide):

1. Without executing the script, explain what happens when you run the script. You should do this by reading the code. As a group, trace execution of the script by hand.
2. Name any places you see context managers or decorators in the script. What do these context managers or decorators do? If part of a library, look up the library documentation and write what the library is used for.
3. This script takes command line arguments (a topic this week!). How does this script get command line arguments from the user?

## Debugging

| Command | Description |
| ------- | ----------- |
| `n`     | Step to the next line of code. |
| `s`     | Step into a function call. |
| `c`     | Continue execution until the next breakpoint. |
| `p`     | Print the value of a variable. |
| `l`     | List the code around the current line. |
| `q`     | Quit the debugger. |

1. First, place a breakpoint at the top of the module and execute `python main.py` in your terminal. What happens? Press the `l` key to list the code around the current line.You will see an output similar to the following:

```
2     Script for retrieving molecule from PubChem and drawing it using networkX
3     """
4  
5     breakpoint()
6  
7  -> import networkx as nx
8     import matplotlib.pyplot as plt
9  
10     # This will be for coloring nodes in the visualization
11     color_map = {"C": "gray", "O": "red", "N": "blue", "H": "#D3D3D3"}
12   
```

In this example, line 7 (`import networkx as nx`) is the line that will be executed next (for example if you press the `n` key). The surrounding lines provide context, making it easier to understand the flow of the code.

You can exit the debugger by typing `q` and pressing enter.

2. What happens if you place a breakpoint under `if __name__ == "__main__"` and execute `python main.py` in your terminal?

Next, try running the command mentioned above for visualizing caffeine. 

```bash
python main.py caffeine
```

This will result in an error. What is the error? Place a breakpoint in the code to examine the variables where this goes wrong. Remember that you can use introspection functions like `dir()` and `type()` to examine the variables in the debugger.

3. Where did you place the breakpoint in your code?
4. What commands did you use after placing the breakpoint to step through program execution? List each command and what you examined.
5. How would your breakpoint placement have changed if the function causing the problem were in an external library? Why?

Use your Python knowledge to fix this problem. There will still be another error in the code that you have to fix after fixing this initial problem. You might find it useful to use the debugger to solve this problem as well. Describe your debugging process in your README.

## Code Improvement
What else can be done to improve this code? You might consider documentation, commenting, and error handling. For example, what will happen if someone enters a molecule name that is not found on PubChem? Add error handling for this case. 

  
