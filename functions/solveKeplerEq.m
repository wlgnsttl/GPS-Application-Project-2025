function Ek = solveKeplerEq(Mk, e)
%케플러 방정식의 항 Mk, e를 입력받아 Ek 도출(Newton–Raphson method) 
    
    x0 = Mk; % 초기값은 M으로 설정

    thr_E = 1e-10; %임계값
    max_iter = 100; %최대 반복횟수

    numb_iter = 0;
    Err = 100;

    while Err > thr_E && numb_iter < max_iter
        numb_iter = numb_iter + 1;

        fx0 = x0 - e*sin(x0) - Mk;
        dfx0 = 1 - e*cos(x0);

        x1 = x0 - fx0/dfx0;

        Err = abs(x1-x0);
        
        x0 = x1;
    end
    Ek = x1;
end
