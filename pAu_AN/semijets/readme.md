# Semi inclusive jets:

To do:

 + [x] Get the cross section for the events:
```c++
    static constexpr double f_n10xPyth6[9] {
        3920155 , 2101168 , 1187058 , 1695141 ,
        4967075 , 1797387 , 260676 , 261926 , 262366
    };
    static constexpr double f_XsecPyth6[9] {
        0.107509,   0.0190967,  0.00475202,
        0.00198812, 0.000361282, 9.65463E-06,
        4.71077E-07, 2.68464E-08, 1.38211E-09
    };
    static double w10xPyth6(int i) {
        return f_XsecPyth6[i] / f_n10xPyth6[i];
    };
```

 + get the lead jet count, and find the turnon effect

 + look at grid: `/tier2/home/groups/rhi/vverkest/pAu_GracePythEmb/src/third.cc` -- will get jets
 + [ ] look at the turn on for full and charged jets w.r.t. the trigger tower


 just check how the sets look and plot


------------------------------------------------------------------------------------------------------------------
### Code instructions
------------------------------------------------------------------------------------------------------------------
user runs:               comments(//), inputs(<-), output (->) 
                         (Note: unless otherwise indicated, `.x tu_load_libs` and `.L loc_libs.h` are assuemd.
---------------------    -----------------------------------------------------------------------------------------
`scp wsu:/tier2/home/groups/rhi/vverkest/pAu_GracePythEmb/out-data/jetemb/hadd.root hadd_jetemb.root`


