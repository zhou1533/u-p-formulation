  clear;
  clc; 

   p_analytical =[ 0, 4.87391e-45, 0.486438, 0.972876, 0.972876, 0.972876, 0.972876, 0.972876, 0.972878, ...
    0.972885, 0.972903, 0.972943, 0.973019, 0.973148, 0.973348, 0.973638, 0.974038, 0.974565,...
    0.975232, 0.976054, 0.977041, 0.978199, 0.979535, 0.981052, 0.98275, 0.984629, 0.986688, ...
    0.988921, 0.991327, 0.993898, 0.996629, 0.999515, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,...
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];


  p_ImprovedEuler=[ 0   0.0041    0.0671    0.2493    0.5073    0.7458    0.9035    0.9786    1.0037    1.0119   ...
      1.0200    1.0295    1.0350    1.0316 1.0181    0.9963    0.9702    0.9439    0.9214    0.9056    0.8981   ...
      0.8996    0.9092    0.9253    0.9452    0.9661    0.9851    0.9999 1.0089    1.0114    1.0077    0.9989  ...
      0.9869    0.9737    0.9615    0.9520    0.9467    0.9462    0.9503    0.9583    0.9689    0.9805 0.9915  ...
      1.0003    1.0059    1.0077    1.0057    1.0004    0.9928    0.9841    0.9757    0.9688    0.9643    0.9627 ...
      0.9642    0.9683  0.9742    0.9811    0.9877    0.9932    0.9966    0.9976    0.9961    0.9922    0.9866  ...
      0.9802    0.9737    0.9682    0.9641    0.9621 0.9622    0.9642    0.9677    0.9719    0.9761    0.9796   ...
      0.9818    0.9822    0.9809    0.9778    0.9734];
 

   p_Euler= [0.0005   -0.0675    0.1713    0.7401    1.0378    1.1069    1.0720    0.9964    0.9368    0.9245    0.9512    0.9845    ...
        0.9959    0.9800 0.9530    0.9363    0.9402    0.9578    0.9735    0.9758    0.9655    0.9529    0.9485    0.9552  ...
        0.9670    0.9747    0.9736    0.966 0.9601    0.9603    0.9668    0.9745    0.9779    0.9754    0.9703    0.9675  ...
        0.9696    0.9749    0.9796    0.9805    0.9777    0.9743 0.9735    0.9761    0.9801    0.9827    0.9822    0.9796 ...
        0.9775    0.9778    0.9803    0.9832    0.9842    0.9830    0.9809    0.9798 0.9806    0.9828    0.9845    0.9846  ...
        0.9832    0.9815    0.9811    0.9822    0.9838    0.9847    0.9842    0.9828    0.9817    0.9817 0.9827    0.9839  ...
        0.9841    0.9833    0.9820    0.9813    0.9816    0.9825    0.9831    0.9829    0.9819];
    
   p_Zienkiewics= [ 0    0.0049    0.0907    0.2277    0.3676    0.4883    0.5821    0.6476    0.6870    0.7058    0.7120    0.7145 ...
        0.7205    0.7340 0.7547    0.7791    0.8018    0.8180    0.8253    0.8239    0.8171    0.8090    0.8038 ...
        0.8035    0.8083    0.8160    0.8236    0.8283 0.8287    0.8249    0.8188    0.8125    0.8083    0.8072 ...
        0.8088    0.8118    0.8144    0.8150    0.8130    0.8087    0.8035    0.7986 0.7954    0.7941    0.7944 ...
        0.7953    0.7956    0.7945    0.7917    0.7876    0.7832    0.7794    0.7768    0.7756    0.7752    0.775 ...
        0.7742    0.7723    0.7692    0.7656    0.7619    0.7588    0.7567    0.7555    0.7548    0.7540    0.7526  ...
        0.7503    0.7474    0.7442 0.7411    0.7386    0.7369    0.7357    0.7347    0.7335    0.7318    0.7295...
        0.7268    0.7240    0.7214];
    
  

       u_analytical= [0  -12.9161  -19.7282  -25.5432  -30.8407  -35.7953  -40.5182  -45.0656  -49.4734  -53.7548  -57.9475  -62.0517 ...
             -66.0930  -70.0706 -73.9850  -77.8488  -81.6619  -85.4628  -89.2000  -92.9119  -96.5859 -100.2476 -103.8709 ...
             -107.4690 -111.0417 -114.5891 -118.1243 -121.6337 -125.1306 -128.6021 -132.0609 -135.5075 -138.9410 -142.3493 ...
             -145.7449 -149.1408 -152.5111 -155.8687 -159.2137 -162.5587 -165.8913 -169.2110 -172.5180 -175.8124 -179.0941 ...
             -182.3888 -185.6579 -188.9269 -192.1707 -195.4274 -198.6712 -201.9023 -205.1207 -208.3782 -211.5579 -214.7377 ...
              -217.9048 -221.2020 -224.3565 -227.5110 -230.6655 -233.9373 -237.0665 -240.1957 -243.3248 -246.4413 -249.6879 ...
              -252.7918 -255.9083 -258.9995 -262.2334 -265.3246 -268.4031 -271.4943 -274.5729 -277.7688 -280.8474 -283.9133 ...
              -286.9792 -290.1625 -293.2157];
      

     u_Euler= [ -0.0147   -7.9199  -20.4657  -30.2272  -35.4554  -37.3718  -39.0418  -42.9364  -49.3012  -56.3734  -62.0196  -65.4224 ...
       -67.5585  -70.2476 -74.6520  -80.4367  -86.1767  -90.5816  -93.5073  -95.9794  -99.2887 -103.9452 -109.3291 -114.2673 ...
       -118.0100 -120.8048 -123.6342 -127.3817 -132.1287 -137.1379 -141.4953 -144.8567 -147.6937 -150.8685 -154.9173 -159.6355 ...
       -164.2819 -168.2060 -171.3694 -174.3464 -177.8292 -182.0619 -186.6651 -190.9726 -194.5854 -197.6824 -200.8491 -204.5844 ...
       -208.8998 -213.3321 -217.3405 -220.7488 -223.8732 -227.2440 -231.1662 -235.4809  -239.7143 -243.4760 -246.7726 -249.9867 ...
       -253.5499 -257.5904 -261.8408 -265.8714 -269.4463 -272.7080 -276.0464 -279.7721 -283.8664 -288.0096 -291.8524 -295.3021 ...
       -298.5878 -302.0623 -305.9097 -310.0020 -314.0147 -317.6990 -321.0812 -324.4317 -328.0356];
    


   u_Zienkiewics=[-0.0146   -4.8217  -15.5807  -28.2323  -39.6761  -48.1078  -53.1302  -55.5506  -56.8900  -58.7494  -62.2384 ...
       -67.6469  -74.4474  -81.5934 -87.9766  -92.8617  -96.1374  -98.3073 -100.2440 -102.8212 -106.5735 -111.5133 ...
       -117.1654 -122.7886 -127.6880 -131.4871 -134.2578 -136.4647 -138.7579 -141.7042 -145.5671 -150.2196 -155.2166 ...
       -159.9904 -164.0841 -167.3284 -169.8951 -172.2114 -174.7767 -177.9593 -181.8557 -186.2651 -190.7841 -194.9780 ...
       -198.5567 -201.4850 -203.9852 -206.4368 -209.2175 -212.5513 -216.4259 -220.6070 -224.7427 -228.5095 -231.7417 -234.4917 ...
-236.9998 -239.5905 -242.5385 -245.9599 -249.7720 -253.7345 -257.5523 -260.9968 -263.9955 -266.6551 -269.2122 -271.9330 ...
-275.0047 -278.4633 -282.1846 -285.9400 -289.4930 -292.6948 -295.5409 -298.1684 -300.7950 -303.6290 -306.7850 -310.2415 -313.8538];


   u_ImprovedEuler=[-0.0085   -3.0209  -10.2918  -19.5380  -29.0731  -37.9973  -45.9064  -52.6183  -58.0598  -62.2628  -65.3824  -67.6913 ...
       -69.5415  -71.3108 -73.3482  -75.9278  -79.2186  -83.2696  -88.0117  -93.2756  -98.8216 -104.3776 -109.6798 -114.5102 ...
       -118.7257 -122.2751 -125.2032 -127.6397 -129.7775 -131.8422 -134.0587 -136.6184 -139.6529 -143.2173 -147.2847 ...
       -151.7532 -156.4634 -161.2233 -165.8380 -170.1377 -174.0023 -177.3776 -180.2812 -182.7981 -185.0662 -187.2551 -189.5400...
       -192.0761 -194.9761 -198.2952 -202.0237 -206.0910 -210.3764 -214.7291 -218.9906 -223.0180 -226.7057 -229.9991 ...
       -232.9024 -235.4761 -237.8274 -240.0934 -242.4197 -244.9395 -247.7537 -250.9168 -254.4297 -258.2409 -262.2547 -266.3464 ...
     -270.3814 -274.2354 -277.8120 -281.0573 -283.9667 -286.5846 -288.9973 -291.3194 -293.6764 -296.1856 -298.9395];

    
    sum=zeros(4,2);
    Goodness_of_Fit=zeros(3,2);
    
        for j=1:81
        
          sum(1,1)=sum(1,1)+(p_ImprovedEuler(j)-p_analytical(j))^2;
          sum(2,1)=sum(2,1)+(p_Euler(j)-p_analytical(j))^2;
          sum(3,1)=sum(3,1)+(p_Zienkiewics(j)-p_analytical(j))^2;
          sum(4,1)=sum(4,1)+p_analytical(j)^2;
          
          
          sum(1,2)=sum(1,2)+(u_ImprovedEuler(j)-u_analytical(j))^2;
          sum(2,2)=sum(2,2)+(u_Euler(j)-u_analytical(j))^2;
          sum(3,2)=sum(3,2)+(u_Zienkiewics(j)-u_analytical(j))^2; 
          sum(4,2)=sum(4,2)+u_analytical(j)^2;
          
          
             
        end
    %������;�������� ���� ���������õ��Ǿ�������� 
    
        for i=1:3
    
              Goodness_of_Fit(i,1)=1-sqrt((sum(i,1))./sum(4,1));%��϶ѹ��
              Goodness_of_Fit(i,2)=1-sqrt((sum(i,2))./sum(4,2));%λ��
    
              
              
        end
       Goodness_of_Fit %˳��ΪImprovedEuler Euler Zienkewiks
        
        
        
        
        
        
        
   
    