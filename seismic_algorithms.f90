module seismic_algorithms
   implicit none
   public

contains

   subroutine semblance(suCMPdata, offsets, velcoer, t0_data, dt, nt, ntraces, nvelcoer, coermatrix)
      integer, intent(in) :: nt, ntraces, nvelcoer
      real, intent(in) :: suCMPdata(nt, ntraces)
      real, intent(in) :: offsets(ntraces), velcoer(nvelcoer)
      real, intent(in) :: t0_data, dt
      real, intent(out) :: coermatrix(nt, nvelcoer)

      real, save :: timewincoer = 0.060
      integer :: i, j, k, iw, nw, iv, it
      real :: t0, v, hs, tshift, ftvt, frac
      real :: amplitude, sumamplitude, sumenergia, energiasomatracos, energiatracos
      real, allocatable :: t(:)

      allocate (t(ntraces))

      ! Calculando os indices para a janela temporal
      nw = nint((timewincoer/dt + 1)/2)  ! Numero de pontos para a metade da janela

      do iv = 1, nvelcoer
         v = velcoer(iv)
         do i = 2, nt - 1
            t0 = t0_data + (i - 1)*dt
            k = 0
            do j = 1, ntraces
               k = k + 1
               hs = offsets(j)/v
               t(k) = sqrt(t0*t0 + hs*hs)
            end do
            energiasomatracos = 0.0
            energiatracos = 0.0
            do iw = -nw, nw, 1
               tshift = dt*iw
               sumamplitude = 0.0
               sumenergia = 0.0
               k = 0
               do j = 1, ntraces
                  k = k + 1
                  ftvt = ((t(k) + tshift) - t0_data)/dt
                  it = 1 + int(ftvt)
                  frac = ftvt - int(ftvt)
                  if (it >= 1 .and. it <= (nt - 1)) then
                     amplitude = (1.-frac)*suCMPdata(it, j) + frac*suCMPdata(it + 1, j)
                     sumamplitude = sumamplitude + amplitude
                     sumenergia = sumenergia + amplitude*amplitude
                  end if
               end do
               energiasomatracos = energiasomatracos + sumamplitude*sumamplitude
               energiatracos = energiatracos + sumenergia
            end do
            coermatrix(i, iv) = energiasomatracos/(ntraces*energiatracos)
         end do
      end do
   end subroutine semblance

   SUBROUTINE read_velpicks(inunit, filename, cdp, tnmo_picks, vnmo_picks, nmo_picks_size)
      !Read from a text file the time-velocity pairs of a CDP that were saved by speed analysis
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: inunit, cdp
      CHARACTER(LEN=80), INTENT(IN)  :: filename
      integer, intent(in) :: nmo_picks_size
      REAL, INTENT(OUT) :: tnmo_picks(nmo_picks_size), vnmo_picks(nmo_picks_size)
      INTEGER :: openstatus, inputstatus, checkstatus
      INTEGER :: i, j, k, tamanho_str, n_separator, n_cdps, icdp_loc, n_tnmo
      CHARACTER(len=160) :: cdp_string, tnmo_string, vnmo_string, tmp_string
      CHARACTER :: separator = ","
      INTEGER, DIMENSION(:), ALLOCATABLE :: cdp_list

      OPEN (UNIT=inunit, FILE=filename, ACTION='READ', IOSTAT=openstatus)
      IF (.NOT. (openstatus .EQ. 0)) STOP "read_velpicks: Velocity picks file open failed!"

      READ (inunit, "(A)", iostat=inputstatus) cdp_string

      tamanho_str = LEN_TRIM(ADJUSTL(cdp_string)); 
      tmp_string = cdp_string(5:tamanho_str); ! remove o texto "cdp=" da lista

      n_separator = 0
      DO i = 1, LEN_TRIM(tmp_string)
         IF (tmp_string(i:i) == separator) THEN
            n_separator = n_separator + 1; 
         END IF
      END DO
      n_cdps = n_separator + 1; 
      ALLOCATE (cdp_list(n_cdps), STAT=checkstatus)
      IF (checkstatus > 0) STOP "read_velpicks: Fail! allocating cdp_list"

      READ (tmp_string, *) cdp_list ! lendo os cdps da lista a partir do string

      ! Localizando o indice do CDP atual na lista de CDPs
      icdp_loc = 0; 
      DO i = 1, n_cdps
         IF (cdp_list(i) == cdp) THEN
            icdp_loc = i; 
            EXIT  ! sai do laco
         END IF
      END DO

      IF (icdp_loc == 0) THEN
         WRITE (*, *) "Error! There is no actual CDP in the file with the velocity picks"
         STOP ! interrompe o programa pq nao encontrou do CDP atual na lista de CDPs dos picks
      END IF

      ! Lendo do arquivo o par de linhas tnmo e vnmo correspondentes ao cdp atual
      DO i = 1, n_cdps
         READ (inunit, "(A)", iostat=inputstatus) tnmo_string
         READ (inunit, "(A)", iostat=inputstatus) vnmo_string
         IF (cdp_list(i) == cdp) EXIT
      END DO

      tamanho_str = LEN_TRIM(ADJUSTL(tnmo_string)); 
      tmp_string = tnmo_string(6:tamanho_str); ! remove o texto "tnmo=" da lista

      n_separator = 0
      DO i = 1, LEN_TRIM(tmp_string)
         IF (tmp_string(i:i) == separator) THEN
            n_separator = n_separator + 1; 
         END IF
      END DO
      n_tnmo = n_separator + 1; 
      ! ALLOCATE (tnmo_picks(n_tnmo), vnmo_picks(n_tnmo), STAT=checkstatus)
      ! IF (checkstatus > 0) STOP "Read_velpicks: Fail! allocating tnmo_picks and vnmo_picks"

      READ (tmp_string, *) tnmo_picks ! lendo os tempos nmo a partir do string

      tamanho_str = LEN_TRIM(ADJUSTL(vnmo_string)); 
      tmp_string = vnmo_string(6:tamanho_str); ! remove o texto "vnmo=" da lista

      READ (tmp_string, *) vnmo_picks ! lendo as velocidades nmo a partir do string

   END SUBROUTINE read_velpicks

   SUBROUTINE velpicks_to_trace(npicks, tnmo, vnmo, t0_data, dt, nt, vnmo_trace)
      !Interpolates the time-velocity pair of the current CDP for each time sample of a CDP trace
      implicit none
      integer, intent(in) :: npicks
      real, intent(in) :: tnmo(npicks), vnmo(npicks)
      real, intent(in) :: t0_data, dt
      integer, intent(in) :: nt
      real, intent(out) :: vnmo_trace(nt)
      integer :: i, j, it0, it1, it2

      ! extrapolando a primeira velocidade dos picks para todas as amostras acima do primeiro tnmo
      it0 = INT((tnmo(1) - t0_data)/dt) + 1; 
      vnmo_trace(1:it0) = vnmo(1); 
      ! Interpolacao linear entre os picks das velocidades primeira, segunda, ate alguma velocidade
      DO j = 1, npicks - 1
         it1 = INT((tnmo(j) - t0_data)/dt) + 1; 
         it2 = INT((tnmo(j + 1) - t0_data)/dt) + 1; 
         DO i = it1 + 1, it2
            vnmo_trace(i) = vnmo(j) + (i - it1)*(vnmo(j + 1) - vnmo(j))/(it2 - it1); 
         END DO
      END DO

      ! extrapolando a ultima velocidade dos picks para todas as amostras abaixo do ultimo tnmo
      it0 = INT((tnmo(npicks) - t0_data)/dt) + 1; 
      vnmo_trace(it0:nt) = vnmo(npicks); 
   END SUBROUTINE velpicks_to_trace

   SUBROUTINE apply_nmo(ntracesCMP, nt, t0_data, dt, CMPdata, offsets, vnmo_trace, smute, CMPdata_nmo)
      implicit none
      integer, intent(in) :: ntracesCMP, nt
      real, intent(in)  :: t0_data, dt
      real, intent(in) :: CMPdata(nt, ntracesCMP)
      real, intent(in) :: offsets(ntracesCMP), vnmo_trace(nt)
      real, intent(inout) :: smute
      real, intent(out) :: CMPdata_nmo(nt, ntracesCMP)
      INTEGER :: i, j, k, it, ih, itvt, it0
      REAL    :: t, t0, hs, ftvt, frac, stretchScale, amplitude, stretchFactor, ft0
      INTEGER, PARAMETER :: lmute = 13 ! Length (in samples) of linear ramp for stretch mute
      REAL    :: dtaper(lmute)

      DO it = 1, nt
         t0 = t0_data + (it - 1)*dt; 
         DO ih = 1, ntracesCMP
            hs = offsets(ih)/vnmo_trace(it); 
            t = SQRT(t0*t0 + hs*hs); 
            ftvt = (t - t0_data)/dt; 
            itvt = 1 + INT(ftvt)
            frac = ftvt - INT(ftvt)
            stretchScale = t0/t; ! To multiply the NMO-corrected samples by NMO stretch factor
            IF (itvt < 1 .OR. itvt > nt) CYCLE
            ! aplica nmo
            amplitude = (1.-frac)*CMPdata(itvt, ih) + frac*CMPdata(itvt + 1, ih); 
            CMPdata_nmo(it, ih) = amplitude*stretchScale; 
         END DO ! loop for offsets
      END DO ! loop for time

      ! Aplicacao de mute automatico usando fator stretch nos tracos com NMO e os eventos estirados
      CALL muteTaperRamp(lmute, dtaper) ! calcula os coeficientes (entre 0 e 1) da funcao taper
      IF (smute < 1.5) smute = 1.5  ! revisar pq da e
      DO it = 1, nt
         t0 = t0_data + (it - 1)*dt; 
         DO ih = 1, ntracesCMP
            hs = offsets(ih)/vnmo_trace(it); 
            t = SQRT(t0*t0 + hs*hs); 
            stretchFactor = (t - t0)/t0 + 1; ! Samples with NMO stretch exceeding smute are zeroed
            IF (stretchFactor >= smute) THEN
               ft0 = (t0 - t0_data)/dt; 
               it0 = 1 + INT(ft0); 
               CMPdata_nmo(1:it0, ih) = CMPdata_nmo(1:it0, ih)*0.; 
               DO k = 1, lmute
                  i = it0 + k; 
                  CMPdata_nmo(i, ih) = CMPdata_nmo(i, ih)*dtaper(k); 
               END DO
            END IF
         END DO ! loop for offsets
      END DO ! loop for time

   END SUBROUTINE apply_nmo

   SUBROUTINE muteTaperRamp(nramp, dtaper)
      ! Calculates the coefficients of the linear taper function for muting
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nramp
      REAL, INTENT(OUT), DIMENSION(nramp) :: dtaper
      INTEGER :: i, j
      j = 0; 
      DO i = nramp, 1, -1
         j = j + 1; 
         dtaper(j) = REAL(nramp - i)/REAL(nramp); 
      END DO
   END SUBROUTINE muteTaperRamp

end module seismic_algorithms
