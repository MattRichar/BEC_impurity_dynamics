SRC= Gen_RK4.cpp Double_Well_functions.cpp Schrodinger.cpp
EXENAME=SchrodingerN50.exe
CMP=c++ -std=c++14

#$at means I run whatever I am doing for the target which is the first statement written. 

OBJ=$(SRC:.cpp=.o) #This command sets OBJ as a variable which makes all of the files listed in SRC as object files (the file that compiles).

.phony: r c

$(EXENAME):$(OBJ) #This term makes the executebale file from the objects. I specified what I wanted the exectuable file to be called, now I need to list all the files I want it to run with. 
	$(CMP) $(OBJ) -o3 -march=native -Wall -o $@ 
#This is a similar format as to how you normally compile a program, compiler   objects    -o flag. The $at indicates what target I want to use this for which in this case is the EXENAME

$(OBJ):%.o:%.cpp #To make anything inside of the object list I need to specify how to make this file. $lessthan means that I use everything from the .c then state -c we are compling, -o linking into whatever is on the right
	$(CMP) $< -c -o  $@ 

r:$(EXENAME)
	./$(EXENAME) 
#This command allows me to type "make r" into linux and it will run the executable file which all the other complies and codes. 
c:
	-rm -f *.o #This command allows me to type "make c" and it will remove all of the object files from the f90s. I can also type "make c r" and it will delete all the object files then recreate them. 
# You do not need to run make again unless you make any changes. 
