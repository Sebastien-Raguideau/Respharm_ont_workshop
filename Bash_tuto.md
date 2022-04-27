


## Basics of bioinformatic : Introduction to command lines
Most workflows in bioinformatics are built with command lines in a terminal. The way to interact with files is different from graphical interface, you can't click on anything and need to type a command for the slightest operation, even changing directory. However, terminal is much more versatile than graphical interface : you dispose of various tools to explore your files, do simple text manipulation or just write simple scripts to automatize a task.   
 

### Navigation using command lines
A command is a script/program/application .... installed on your device. To use them, you need to type their name followed by argument and options. Example : 

    Command -option argument
All command presented below are always present on bash terminals.


**Get information about a function** : *man*

**List all files and folders** : *ls*

Try to use **man**  on **ls** A lot of different options are available. Try for yourself  the following : 
 - ls -l
 - ls -lh
 - ls -a
 - ls -alh

Which one are folders, which one are files? 

**Change directory** : *cd*

Try to go into  Dummy :

    cd Dummy
It's not working for some reason? Please take the habit of reading errors message they are often informative.

Go to directory data :

    cd data
You should now see a path on your terminal :

*ubuntu@Respharm:~/data*

This indicate where you are : 

 - **~** is your home folder
 - **data** is the current folder you are in

Take the habit of looking at this to keep track of which folder you are in. 

You can go back to the parent directory by using  :

    cd   ../

You can also use a full such as  :

    cd data/mystery_16S
and you can go back more than once 

    cd ../../

Finally if you are ever lost just type `cd` and you'll be back to your home folder

### create folder and global variable
Bash is not just a way to replace the graphic interface and look at file. It is also a language in which you can for instance define variables:
```bash
export DATA=/home/ubuntu/data
```
We just assigned a value to DATA which can now be referred to. Try the following what happens
```bash
    echo $DATA
    ls $DATA
    cd $DATA
```
 
 ### Working with multiple person on the same folder

As the number of VM is limited we will need you to all work in separate folders. To simplify the rest of the tutorial please create a folder with your name using the command mkdir:

obviously replace \<yourname\> by ... your name 
```bash
cd 
mkdir <yourname>
```  
Then we will store the path to your personal directory inside a variable as before:
```bash
export DATA=/home/ubuntu/<yourname>
```
One issue is that these are only lasting for the current session.

  ### Simple files exploration 
 **Look into a file**  
```bash
cd  $DATA/mystery_16S
```

How big are the files in this folder ? What would happen if you were to try and open then in a text editor?
Look into them by using 
 - less 
 - more 
 - head
 - tail
 
Type q to  quit.

**Number of line in a file** 
*wc -l file*

    wc -l  BC01.fastq.fq

 A fastq entry contain 4 lines for each read. How many reads are there?  
 
**Search for a pattern**   
*grep pattern file*
**grep** will return any line of the file containing the pattern. 
Lets try to see if we can find some absurds reads :

    grep AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA BC01.fastq.fq

Which reads are they? Lets use the option -B 1 to show the line just before the pattern.

    grep -B 1 AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA BC01.fastq.fq

### Combination of commands
You can combine commands by sending the output of a command line as input to a second one using the symbol pipe **|**.
Lets try to count the number of stop codon in a Fastq file. 

    grep TAG BC01.fastq.fq | wc -l


### Text edition 

**Writing**  
*Command line > file*

Using the symbol **>** will allow to write the output of previous command into a file. 

     grep TAG BC01.fastq.fq > Codons_TAG.txt 

**Text editors** 

A few text editors are available on the terminal, 

 - vim
 - emacs
 - nano

Legend states that after spending a few months learning how to use them , vim and emacs are fantastics text editors!
For everyone else, nano exist.

    nano Codons_TAG.txt



**Remove file**  
*rm file*
Be careful with this one.

    rm Codons_TAG.txt

### Writing a simple script 
The terminal allow you to write scripts. You can for instance loop through files and  apply the same treatment. 
Let's write a simple script which will loop through all folders of Data, and count the number of erroneous reads. It is going to take time executing, so feel free to interrupt it at any time using key combination ctrl + c
```Bash
    cd $DATA
    for file in */*/*
    do 
    	echo $file
    	echo $file >> Number_of_AAA_in_all_Data_files.txt
    	grep AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA $file | wc -l >> Number_of_AAA_in_all_Data_files.txt 
    done
```
 - **\*** is called a wild-card : it can replace anything, here it is used to select everything inside all folders insides all folders. 
 - **file** is variable, to access the value of the variable you need to add **$**
 - **echo** is a command used to print a variable, so it is responsible for printing the value of **file** on screen, and in Number_of_AAA_in_all_Data_files.txt.
 - We use **>>** instead of **>**, it allows to append at the end of a file, while **>**  would recreate a file each time
 

