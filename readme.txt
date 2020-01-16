Prince, Keith

As I’d hoped here is some python code which will produce the real-time wave analysis

I’ve written it in a fairly modular way so the structure of the main code is

1)	Defining some parameters, you shouldn’t need to change any of these, but in time you may want to make some of your
geophysical constants consistent with what you use in the model, it shouldn’t make an significant difference to the
output though

2)	Read in 90days of data @6 hour resolution, here I’ve done it on a one degree grid, but I think the code should be
able to work out the resolution, there is no need for higher resolution though as we’re throwing away anything which
is shorter than wave No 40 anyway.

We are expecting (and what we have tested it with is) the last 29 points to be the analysis time (AT) and the next
7 days of forecast such that the last point is T+168, and the preceding days run from (AT-1986:At-6)

The data is dimension (in python space) as [time,pressure,lat,lon] it’s coded to work on more than one pressure level,
but we’ve only really been looking at the low-level (850hPa) as that’s what we’ve used for all our other analysis.
It will work for one pressure level providing the pressure coordinate is include e.g. n_levels =1

The code only needs latitudes to +/- 24 with an equatorial point and the read_data routine as written reads global
fields and extracts that particular subset and should be agnostic to the ordering of latitudes (N-S or S-N)

You may want to replace this with a routine to get the data from the analysis and forecasts or write a script that
will produce a file in the same format.

Now the actual projection onto the waves begins
3)	Convert u and z to new variables: these are the appropriate variables for projecting on to wave modes following
the theory
4)	Fourier transform in time and space: Standard numpy fft routines
5)	Project onto the wave modes: The gory stuff is all in here, it could be extended and/or speeded up, but I would
advise leaving it as it is at the moment, this is place where all the bugs were (even just 2 hours ago), takes the
spectral q,r, and v as inputs and returns spectral u,v,z for each wave mode.
Note as it stands qf,rf, and vf are re-ordered in this routine, but I think as they’re not returned they should
come out unchanged (check if you want to make more use of them)
6)	Inverse FFT: standard fft routines again, returns the physical u,v,z for each wave

Currently this is all coded to do Kelvin Waves; Westward Moving Mixed Rossby-Gravity Waves;
and n=1,2 Rossby Waves (R1,R2)

You can reduce the number of waves you output, but it want substantially reduce the amount of run time as most of the
calculation is unaware of that choice (but could in principle be made to be)

7)	Write out the data: The wave data is all written to one file, as five dimensional fields for u,v,z in python space
the dimensions are [wave_type,time,pressure,lat,lon) the ordering of the wave dimensions is as defined in the code
(currently, [‘Kelvin’,’WMRG’,’R1’,’R2’])
Currently we’re writing out the whole 90 days of data, but realistically you will probably only want to write out say
(T-4:T+7) so that you can plot an animation of the last few days and the forecast.
Note that the Kelvin Wave v is 0 by definition

Note: Before plotting the data on the analysis you should do the bias correction step described in Yang’s report which
removes the mean of the last 30 days lead time dependent analysis; that is
a.	Produce a wave analysis for 30 consecutive days (e.g. for forecast starting on the 1st of October you would also
need to do the analysis for the whole of September)
b.	For each time remove the lead time dependent climatology of the previous 30 days ; i.e. average the 30 T+6 waves
and plot the forecast T+6 wave as an anomaly from that (and so on)
In operational practice this will mean keeping the wave analysis data so that the bias correction can be done, and
hence I would suggest only keeping e.g. T-4:T+7 to reduce the data volume

When plotting these waves on a map we have used the convergence to plot the Kelvin Waves and WMRG waves and the
Vorticity for the Rossby Waves (that’s what the plots in the documents you have will have as contours); Sam can
probably let you know what sensible c.i. are for each wave mode to pickout high amplitude waves)

We have also been developing an additional local wave phase diagram similar to WH diagrams for the MJO, and Sam has
recently calculated these for past data, and we can probably work on some code that would take the output of the
wave analysis code and plot them.

I’ve also put a copy of the input and output files (as well as the code) on the ceda website for the project
in the realtime directory @

http://gws-access.ceda.ac.uk/public/ncas_climate/seasia_waves/realtime/

so that you can check things if you want to.

Hope everything  works OK, let me know if you have any questions.

Steve

----------------------------------
Professor Steve Woolnough
National Centre for Atmospheric Science

Department of Meteorology, University of Reading
Earley Gate, PO Box 243, Reading, RG6 6BB, UK
