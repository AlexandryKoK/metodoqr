program metodoqr

!Universidade de Brasilia - Metodos Computacionais A
!Professor: Luiz A. Ribeiro Junior - Instituto de Física
!Aluno: Alexandry Moreira Alves Pinto - 17/0078761
!Problema 3 - Encontrar autovalores utilizando o Metodo QR.
!Github url: 


double precision, dimension(100,100) :: matrix
double precision, dimension(100,100) :: matrixinicial
double precision, dimension(100,100) :: matrixlinha
double precision, dimension(100,100) :: matrixq
double precision, dimension(100,100) :: matrixqt
double precision, dimension(100,100) :: matrixqk
double precision, dimension(100,100) :: matrixqsoma
double precision, dimension(100,100) :: matrixr
real*8 :: modulo, soma, somadois, somatres, novoa
integer :: i,j,k,n,v,p,l,z



write(*,*) "Bem vindo ao programa de autovalores usando o metodo QR."
write(*,*)
write(*,*) "O programa irá ler o arquivo 'entradasqr.txt' "
write(*,*) 
write(*,*) "Digite a ordem da sua matriz (numero inteiro):"
write(*,*)
read (*,*) n
write(*,*)
write(*,*)

open (1, file = 'entradasqr.txt', status = 'old')
do i = 1, n
	do j = 1,n
		read(1,*) matrix(i,j)
		matrixinicial(i,j) = matrix(i,j)
	end do
end do
close(2)



do z=1,50

!Limpando as matrizes para sempre que mudar z, ter o programa inicial
do i=1,n
do j=1,n
matrixq(i,j) = 0
matrixr(i,j) = 0
matrixlinha(i,j) = 0
matrixqk(i,j) = 0
matrixqt(i,j) = 0
matrixqsoma(i,j) = 0
end do
end do

!Para i=1, temos que indicar a coluna 1 da matriz q
do i=1,n
	modulo=0
		do j=1,n
			modulo = modulo + matrix(j,1)**2
		end do
	matrixq(i,1) = matrix(i,1) / sqrt(modulo)
end do

!Indicar que matriz linha 1 = coluna 1 matriz q
do i=1,n
	matrixlinha(i,1) = matrixq(i,1) 
end do


!Dar o primeiro valor para a matriz r
matrixr(1,1) = sqrt(modulo)
somatres=0
do i=1,n
	somatres = somatres + matrixq(i,1) * matrix(i,2)
end do
matrixr(1,2) = somatres


!Para i=2 em diante temos o seguinte algortimo: [a'_i = a_i - Sum(k=1,k=i-1)(q_(k)^(t)a_i)q_k]
!Onde as matrizes a_i são matrizes colunas de A e matrizq São matrizes colunas onde cada elemento
!é qk = a_k' / Modulo(a_k')
do i=2,n

	do k=1,i-1
		soma=0
		do v=1,n
			soma = soma + matrixq(v,k) * matrix(v,i)
		end do
		do v=1,n
			matrixqk(v,k) = soma * matrixq(v,k)
		end do
	end do

	do v=1,n
		do k=1,i-1
			matrixqsoma(v,i) =  matrixqsoma(v,i) + matrixqk(v,k)
		end do
		matrixlinha(v,i) = matrix(v,i) - matrixqsoma(v,i)
	end do

	modulo=0
	do j=1,n
		modulo = modulo + matrixlinha(j,i)**2
	end do
		matrixr(i,i) = sqrt(modulo)
	do j=1,n
		matrixq(j,i) = matrixlinha(j,i) / sqrt(modulo)
	end do



	do k=1,n
	somatres=0
		do p=1,n
			if (k/=i) then
				somatres = somatres + matrixq(p,k) * matrix(p,i)
			end if
		end do
			if (k/=i .and. k<i) then
				matrixr(k,i) = somatres
			end if
	end do


end do



do i=1,n
	do j=1,n
		novoa = 0
			do k=1,n
				novoa = novoa + matrixr(i,k)*matrixq(k,j)
			end do
		matrix(i,j) = novoa
	end do
end do

end do

do i=1,n
	do j=1,n
		if (i==j) then
			write(*,*) "Autovalor λ",i,"=",matrix(i,j)
		end if
	end do
end do

write(*,*)



open(2, file = 'saidaqr.txt', status = 'replace')
write(2,*) "Autovalores aproximados:"
write(2,*)
do i=1,n
	do j=1,n
		if (i==j) then
			write(2,*) "Autovalor λ",i,"=",matrix(i,j)
		end if
	end do
end do

close(2)
write(*,*)
write(*,*) "Salvo em 'saidaqr.txt'"



end program Metodoqr