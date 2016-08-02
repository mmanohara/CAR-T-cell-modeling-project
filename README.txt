README:

This folder contains all of the programs that I used in this project. This doesn’t include the programs where I changed the parameters, but it isn’t hard to do. Instructions on how to use these are below.

These codes are to be run on python. The modules matplotlib and pylab are necessary.

1. Finding the initial condition equation.
	Use the Find constant program.
	The initial condition is assumed to be a Gaussian distribution with major and minor axes representing the dimensions of the region. Each region has the equation Ce^(-(x-xcenter)^2/(xradius)^2 - (y - ycenter)^2/(yradius)^2, where C is a constant that we need to find. First, input the population density data for each cancer species and assign it to a variable RcX, where R is the region letter and X is the cancer species. For example, Ec2 is the density of cancer species 2 in region E. Then, find a carrying capacity k. Divide all of your cancer species by k. Next, measure the dimensions of the regions using an image (this has to be done manually. Look below for more specific instructions.) Input your leftmost x as a list x1, and your rightmost x as a list x2, your topmost y as list y1, and bottommost y as list y2 (measure from top to bottom, not vice versa. Makes life easier). For each region, make a list of the population densities for the cancer from that region (e.g.: A = [Ac1, Ac2, Ac3]). Then, to find the constant for each cancer species in each region, change the region number and the region letter for each region, and run the program. Your final equation will look like:

c1 = C1e^(-(x-xcenter[0])^2/(xradius[0])^2 - (y - ycenter[0])^2/(yradius[0])^2 + C2e^(-(x-xcenter[1])^2/(xradius[1])^2 - (y - ycenter[1])^2/(yradius[1])^2…Cre^(-(x-xcenter[r])^2/(xradius[r])^2 - (y - ycenter[r])^2/(yradius[r])^2

where r is the region number.

The Find constant program contains the data from PBT027. 

A. Measuring the dimensions.
Download pixelstick. Zoom out until the image is completely visible. Measure the length of the long segment (to be denoted L). Measure the length of the short segment (to be denoted S). Calculate L-S/2 (to be denoted B). This value will be added to each of the x values in order to ensure a square space. Using pixel stick, measure the distance from the left edge to the leftmost area of the region, and the same to the rightmost area. Add B to these coordinates. This is x1 and x2. Do this for each region. Next, measure the distance from the upper edge to the upmost area of each region, and the distance to the bottommost area of the region. This is y1 and y2. Afterwards, multiply each by 1/ L to normalize the space to 1x1. Round to the hundredth. Use these in the Find constant program.

2. Using the “find totals” program

You will need to come up with your initial condition equation.	Add the lists of x1, x2, y1, and y2 like in the previous program. Then, replace the initial condition equation with the new initial condition equation. Run the program and you will find the total initial population for each species.

3. Using the “show initial condition” program

Once again, you need to include the x1, x2, y1, and y2 lists. Replace the initial condition equation with the new one. As of now it only includes c1, c2, and c3, but adding cco shouldn’t be too difficult. This program prints the population sizes of the initial condition in addition to showing a distribution of the initial condition.

4. Running simulations

6 programs are there for 6 different treatments. You can adjust the parameters as you wish at the beginning. The program samelocbiCAR contains the most information about variables, so look there for information on the variables. (Also I’m lazy and don’t want to rewrite the comments.) To run it, adjust the parameters as you wish, then replace the initial condition. You will need to actually calculate xcenter, xradius, ycenter, and yradius and put them into the equations, but if you don’t want to you can copy the parts of the “find totals” program which contains the lists and the procedure for generating lists of xcenter, xradius, ycenter, and yradius. Global variables are used so that the number of arguments are decreased. When you run the program, it will generate pics in the same folder as the program, so to isolate the pictures you will have to place the program in a separate folder. 

The 6 treatments all include CAR T-cell doses after every 7 time units. However, the treatments are slightly different for each. Read more in the introduction for each program.

If you need information on any of the programs, my email is mohith_manohara@yahoo.com
