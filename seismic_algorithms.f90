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

end module seismic_algorithms
