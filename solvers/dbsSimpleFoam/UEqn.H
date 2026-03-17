    // Momentum predictor

    tmp<fvVectorMatrix> tUEqn
    (
        1./(eps+SMALL)*fvm::div(phi, U)
      - fvm::laplacian(laminarTransport.nu()/(eps+SMALL),U)
      + fvm::Sp(laminarTransport.nu()*invK,U)
      + fvm::Sp(laminarTransport.nu()*HeleShawCorrection,U)
     ==
        fvModels.source(U)
    );
    fvVectorMatrix& UEqn = tUEqn.ref();


    UEqn.relax();

    fvConstraints.constrain(UEqn);


    if (simple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));

        fvConstraints.constrain(U);
    }
