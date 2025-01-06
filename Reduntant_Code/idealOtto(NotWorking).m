% First average crank angle and take only the power rotation into account
Camean = mean(Ca,2);
idx_half = find(Camean==0);
CaPower = Ca(idx_half:end,:);
VPower = V(idx_half:end,:);
pPower = p(idx_half:end,:);

% Take the mean Ca
CaPower = mean(CaPower);

% average volume over all 69 cycles
Vmean = mean(VPower,2);
[Vmax,idxVmax] = max(Vmean);
[Vmin,idxVmin] = min(Vmean);

% same for pressure
pmean = mean(pPower,2);
pmax = pPower(idxVmin);
pmin = pPower(idxVmax);

Tamb = 293;

R = pmin * Vmax / Tamb;

p1 = pmin;
V1 = Vmax;
T1 = Tamb;

p2 = pmax;
V2 = Vmin;
T2 = p2 * V2 / R;

p4 = pmin;
V4 = Vmax;
T4 = p4 * V4 / R;

V3 = Vmin;
p3 = p4 * V4^1.3 / V3^1.4;
T3 = V3 * p3 / R;

State1 = [V1 p1];
State2 = [V2 p2];
State3 = [V3 p3];
State4 = [V4 p4];

Vscat = [V1 V2 V3 V4];
pscat = [p1 p2 p3 p4];

scatter(Vscat,pscat)