
on 14/02/2021 


### info about potential (important for boundness checking)

###---

3 dimensional case (Q1-Q10-Q29)	==>

# grid Q1 = {-4.,4.}(81 points) 
	   Q10 = {-4.,8.}(121 points) 
	   Q29 = {-4.,4.}(81 points)

	No of energies < -1580. cm^-1 ==> none.
	E_min = -1295.579139 ; location = 161885 (0-based indexing)

	No of points with energy < -1200. ==> 558
	No of points with energy < -1000. ==> 2962

############################################

on 15/02/2021

4 dimensional cases ----

##---

Q1-Q5-Q10-Q29 ==>

(Checking has been done with 11 points for each coordinate starting from -4.,
i.e. 11*11*11*11 grid) ==> (Potential is unbound)

# grid Q1 = {-4.,4.}
# grid Q5 = {-4.,4.}
# grid Q10 = {-4.,8.}
# grid Q29 = {-4.,4.}


# grid Q1 = {-3.8,}
# grid Q5 = {-4.,}
# grid Q10 = {-3.8,}
# grid Q29 = {-3.8,}

No of energies < -1580. cm^-1 ==> none
E_min = 94019.253432 ; location = 14640 (0-based indexing)


# grid Q1 = {-4.,}
# grid Q5 = {-4.,}
# grid Q10 = {-3.8,}
# grid Q29 = {-3.8,}

## Unbound .

# grid Q1 = {-4.,}
# grid Q5 = {-4.}
# grid Q10 = {-3.8}
# grid Q29 = {-3.6}

No of energies < -1580 cm^-1 ==> 0
E_min = 97103.644306 ; location = 14640 (0-based indexing)


##------

Q1-Q7-Q10-Q29 ==>

(11*11*11*11 grid)

# grid Q1 = {-4.}
# grid Q7 = {-4.}
# grid Q10 = {-3.8}
# grid Q29 = {-3.6}

No of energies < -1580 cm^-1 ==> 0
E_min = 93141.886963 ; location = 14640 (0-based indexing)


##------

Q1-Q10-Q13-Q29 ==>

(11*11*11*11 grid)
 
grid Q1 = {-4.}
grid Q10 = {-3.8}
grid Q13 = {-4.}
grid Q29 = {-3.6}

No of energies < -1580 cm^-1 ==> 0
E_min = 92613.140937 ; location = 14640(0-based indexing)

##------

Q1-Q10-Q27-Q29 ==>

(11*11*11*11 grid)

grid Q1 = {-4.}
grid Q10 = {-3.8}
grid Q27 = {-4.}
grid Q29 = {-3.6}

No of energies < -1580 cm^-1 ==> 0
E_min = 113769.481231 ; location = 14640 (0-based indexing)




