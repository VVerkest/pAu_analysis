root -l << EOF
    .L ../loc_lib/noiDict.h

    noiDict d {{1,3,1.24,5.23,"this","my_name","it"}};
    double i = d("1.24");
    cout << i << endl;
    cout << " had 1? " << d["1"] << endl;
    cout << " had 2? " << d["2"] << endl;

    double a = d("3",2);
    cout << " a: " << a << endl;
    double b = d["it"];

    TH1D hg { d("this"),"a;b;c",(int)d("1"),d("1",0),d("1",2) };
    cout << " name " << hg.GetName() << " xaxis bins " << hg.GetXaxis()->GetNbins() << " last bin " << hg.GetXaxis()->GetBinUpEdge(3) << endl;

    cout << d << endl;
    noiDict e {{"samantha","is","a","beautiful","name"}};
    cout << d << endl;
    d += e;
    cout << d << endl;

    noiDict f {{0,0}};
    noiDict g {{1,1}};
    f+= g;
    cout << f << endl;
EOF
