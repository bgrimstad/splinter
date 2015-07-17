
//bool serializePSpline1()
//{
//    // Create new DataTable to manage samples
//    DataTable samples;
//    // Sample function
//    auto x0_vec = linspace(0, 2, 20);
//    auto x1_vec = linspace(0, 2, 20);
//    DenseVector x(2);
//    double y;
//    for (auto x0 : x0_vec)
//    {
//        for (auto x1 : x1_vec)
//        {
//            // Sample function at x
//            x(0) = x0;
//            x(1) = x1;
//            y = sixHumpCamelBack(x);
//            // Store sample
//            samples.addSample(x,y);
//        }
//    }
//    // Build P-spline
//    PSpline pspline(samples);
//    pspline.save("saveTest1.pspline");
//    PSpline loadedPspline("saveTest1.pspline");
//    remove("saveTest1.pspline");
//
//    return compareBSplines(pspline, loadedPspline)
//           && pspline.getLambda() == loadedPspline.getLambda();
//}
//
