# ============================================================
# Fonction : uzawa_solve
# Description : Implementation de l'algorithme d'Uzawa pour
#               resoudre un probleme de point selle issu de
#               l'optimisation sous contraintes d'egalite
#               affines.
#
# Probleme resolu :
#   minimiser  J(u) = 0.5 * u'Au - b'u
#   sous la contrainte  C u = f
#
# Algorithme (iteration k) :
#   1) Etape primale : resoudre A u = b - C' lambda
#   2) Etape duale   : lambda <- lambda + rho * (C u - f)
#
# Condition de convergence sur le pas :
#   0 < rho < 2 * alpha / ||C||^2
#   ou alpha est la constante d'ellipticite de A.
# ============================================================

#' Algorithme d'Uzawa pour contraintes d'egalite affines
#'
#' Resout le probleme de point selle associe a la minimisation
#' de J(u) = 0.5 * u'Au - b'u sous la contrainte C u = f.
#'
#' @param A Matrice symetrique definie positive (n x n).
#' @param b Vecteur du second membre (longueur n).
#' @param C Matrice de contraintes (m x n).
#' @param f Vecteur de contraintes (longueur m).
#' @param rho Pas de mise a jour du multiplicateur.
#'   Doit verifier 0 < rho < 2 * alpha / ||C||^2.
#' @param tol Tolerance sur le residu dual (defaut 1e-6).
#' @param max_iter Nombre maximum d'iterations (defaut 1000).
#'
#' @return Une liste contenant :
#'   \itemize{
#'     \item u : solution primale approchee
#'     \item lambda : multiplicateur de Lagrange approche
#'     \item iterations : nombre d'iterations effectuees
#'     \item converged : TRUE si le critere d'arret est atteint
#'   }
#'
#' @examples
#' A <- diag(2)
#' b <- c(1, 1)
#' C <- matrix(c(1, 1), nrow = 1)
#' f <- 1
#' res <- uzawa_solve(A, b, C, f, rho = 0.5)
#' res$u
#' res$lambda
#'
#' @export
uzawa_solve <- function(A, b, C, f, rho, tol = 1e-6, max_iter = 1000) {

  # --- Verification basique des dimensions ---
  n <- ncol(A)
  m <- nrow(C)

  if (length(b) != n) {
    stop("La longueur de b doit etre egale au nombre de colonnes de A.")
  }
  if (ncol(C) != n) {
    stop("Le nombre de colonnes de C doit etre egal a n.")
  }
  if (length(f) != m) {
    stop("La longueur de f doit etre egale au nombre de lignes de C.")
  }
  if (rho <= 0) {
    stop("Le pas rho doit etre strictement positif.")
  }

  # --- Initialisation du multiplicateur ---
  # On demarre avec lambda = 0 (choix classique).
  lambda <- rep(0, m)

  # Variable pour stocker la solution primale courante.
  u <- rep(0, n)

  # Indicateur de convergence.
  converged <- FALSE

  # Boucle principale de l'algorithme d'Uzawa.
  for (k in seq_len(max_iter)) {

    # --- Etape 1 : resolution primale ---
    # Pour lambda fixe, on minimise le Lagrangien en u :
    #   A u = b - C' lambda
    # On utilise solve() pour resoudre le systeme lineaire.
    rhs <- b - t(C) %*% lambda
    u <- solve(A, rhs)

    # --- Etape 2 : mise a jour duale ---
    # lambda <- lambda + rho * (C u - f)
    # Le residu (C u - f) mesure la violation de contrainte.
    residual <- as.vector(C %*% u - f)
    lambda_new <- lambda + rho * residual

    # --- Critere d'arret ---
    # On compare la variation du multiplicateur a la tolerance.
    delta <- sqrt(sum((lambda_new - lambda)^2))
    if (delta < tol) {
      converged <- TRUE
      lambda <- lambda_new
      break
    }

    # Mise a jour pour l'iteration suivante.
    lambda <- lambda_new
  }

  # --- Avertissement si la convergence n'est pas atteinte ---
  if (!converged) {
    warning("Nombre maximum d'iterations atteint sans convergence.")
  }

  # --- Retour des resultats ---
  list(
    u          = as.vector(u),
    lambda     = as.vector(lambda),
    iterations = k,
    converged  = converged
  )
}
