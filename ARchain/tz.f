      
      CLIGHT = 18000.0
      ECC = 0.308
      SEMI = 7.0D-07
      ECC = 0.99
      SEMI = 2.0D-05
      RAU = 1.0*2.0D+05
      SMU = 6.1D+04
      TAUGR = 1.3D+18*RAU**4/SMU**3
      WRITE (6,1) TAUGR
    1 FORMAT (' TAUGR   ',1P,E10.2)
          ECC2 = ECC**2
          FE = 1.0 + (73.0/24.0 + 37.0*ECC2/96.0)*ECC2
          GE = (1.0 - ECC2)**3.5/FE
          ZX = 3.0D-04
          RATIO = 1.0
*       Replace physical time-scale by N-body units (cf. Lee 1993).
*         TZ = TAUGR*GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZX**3)
          TZ = GE*SEMI**4/(RATIO*(1.0 + RATIO)*ZX**3)
      WRITE (6,3)  SEMI, ZX, TZ
    3 FORMAT (' SEMI ZX TZ  ',1P,3E10.2)
          TZ = 5.0/64.0*CLIGHT**5*TZ
      WRITE (6,6)  TZ
    6 FORMAT (' TZ  ',1P,E10.2)
      STOP
      END 
