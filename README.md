# Analiza pogrešaka RMF modela konačne jezgre

Napravljena je adaptacija RMF procedure za sferične jezgre (Ghambir+90) za primjenu u računu Fisherove metrike koristeći algoritamsku diferencijaciju paketa autograd. 
Konstruirana je Fisherova metrika za set parno-parnih jezgara <sup>8</sup>Be, <sup>12</sup>C, <sup>16</sup>O, <sup>20</sup>Ne, <sup>24</sup>Mg, <sup>28</sup>Si, <sup>32</sup>S, <sup>36</sup>Ar i <sup>40</sup>Ca te je proveden račun pogrešaka koristeći njihove energije vezanja i nabojne radijuse izračunate funkcionalom DD-PC1. Izabrano je da pogreške svakog mjerenja iznose $10\%$ dobivene energije vezanja i nabojnog radijusa.

## Opis procedure
Osnovni kod se bazira na iterativnoj proceduri za određivanje RMF energije vezanja u sfernom sustavu koristeći bazu harmoničkog oscilatora

$R_{n,l}=N_{n,l}L(\xi^2,n,l+0.5)\xi^le^{-\xi^2/2},$

gdje je $\xi$ radijus reskaliran s parametrom $b_0=\sqrt{1.011 A^{1/3}}$, a $L(\xi^2,n,l+0.5)$ pridružen Laguerrov polinom.

Kod rješava Diracovu jednadžbu $H\psi = (\epsilon+m)\psi$ koristeći pojednostavljenu sfernu bazu uz ansatz $\psi=(f(r),ig(r))$, posebno za svaki $(j,\pi)$ blok


Koristeći bazu harmoničkog oscilatora, problem se svodi na blok matricu ako se koristi razvoj 
$f=\sum\limits_{n}^{n_{max}} f_n R_{n,l}$
$g=\sum\limits_{\tilde n}^{\tilde n_{max}} g_{\tilde n} R_{\tilde n,\tilde l}$

gdje su $n_{max}=(N_F-l)/2$, $\tilde n_{max}= N_F+1$,  $l(j,\pi)=l+\pi/2$, $\tilde l(j,\pi)=l-\pi/2$ i $\kappa = \pi (j+1/2)$ . U ovoj bazi problem postaje blok-dijagonalna matrica.

Iz dobivenih svojstvenih vektora, kod računa valne funkcije, s obzirom na sparivanje koristeći sva stanja, svako stanje ulazi u gustoće pomnoženo s faktorima $\nu^2_i$, gdje su faktori $v_i^2$ (s obzirom na degeneraciju)  određeni kemijskim potencijalom $\lambda$ i procjepom $\Delta$:

$v_i^2=(2j_i+1) \frac{1}{2}\left(1-\frac{\epsilon-\lambda}{\sqrt{(\epsilon-\lambda)^2+\Delta^2}}\right).$

Kod pronalazi optimalni kemijski potencijal tako da je zadovoljena jednadžba
$\sum\limits_i v_i^2(\lambda) = N\,ili\, Z.$

Kod tada računa gustoće $\rho_s$, $\rho_v$, $\rho_{tv}$, $\rho_c$, $\Delta\rho_s$, računa DD-PC1 potencijal, uključivo s kulonskim doprinosom, koje koristi u sljedećoj iteraciji te na kraju računa energiju vezanja i nabojni radijus.

## Detalji implementacije

Kako bi se kod mogao povezati s autograd paketom, u autogradu je napisana implementacija Laguerrovih polinoma, strukture blok matrice, trapezna metoda integracije, te funkcije eigh za dijagonalizaciju.
Pri izračunu kulonskog potencijala bilo je potrebno izbjeći singularitet, što se postiglo malom korekcijom nazivnika

$V_c(r)=2\pi e^2\int \frac{\rho_c r'^2 dr' \sin\theta d\theta}{\sqrt{0.01+r^2+r'^2-2rr'\cos\theta}}.$

Kako bi se optimiziralo vrijeme izvođenja, kemijski potencijal je izvrijednjen na fiksnoj mreži, na kojoj se tada tražio minimum razlike $\left|\sum_i v_i^2(\lambda)-(N,Z)\right|$.

Odabran je dovoljno velik broj mjerenja kako bi Fisherova metrika mogla biti nesignularna zbog nedefiniranosti problema, odnosno numerički imala pozitivne svojstvene vrijednosti. Korištenje dviju varijabli (energije vezanja i nabojnog radijusa) je prepolovilo broj jezgara koje je nužno promatrati. Odabrane su parno-parne jezgre kako bi mogli staviti $\Delta\sim 0$, no kod funkcionira za proizvoljnu jezgru.
Zbog specifičnosti implementacije eigh funkcije, jakobijan energije vezanja je morao biti razdvojen na doprinos od promjene konstanti DD-PC1 funkcionala i doprinos od promjene gustoće, koje moraju biti sumirane pri izračunu ukupnog jakobijana po energiji vezanja.

Korišen je DD-PC1 model, reparametriziran fiksnom best-fitting vrijednošću te su unaprijed eliminirani parametri koji su u DD-PC1 modelu egzaktno nula (jer tada nema ovisnosti o $p_\mu$.

Fisherova metrika se računa kao suma po jezgrama, $a$. 

$g_{\mu\nu}=\sum\limits_a \frac{1}{\sigma_{E_a}^2}\partial_\mu E_a \partial_\nu E_a +\sum\limits_a \frac{1}{\sigma_{R_a}^2}\partial_\mu R_a \partial_\nu R_a$
## Stabilnost svojstvenih vrijednosti promjenom broja iteracija

Kod je testiran s obzirom na broj iteracija, počevši od Woods-Saxon potencijala kao početnog uvjeta. Zbog značajnih zahtjeva na radnu memoriju paketa autograd, za izračun jakobijana je korišten mali broj oscilatorskih ljusaka $N_F=3$ te do 15 iteracija. Za veći broj ljusaka, kod je realistično vrtiti s još manjim brojem iteracija. Sami izračun energije vezanja bez poziva jakobijana nema ograničenja na broj ljusaka.
