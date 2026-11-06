# ============================================================
# Algorithme d'Uzawa general
#
# Probleme resolu :
#   min f(x)
#   sous g_i(x) <= 0   (inegalites, i = 1..m)
#   et   h_j(x) = 0    (egalites,   j = 1..p)
#
# f, g, h peuvent etre lineaires ou non lineaires.
# Seule la convexite est supposee pour la convergence.
#
# Principe :
#   Etape primale : minimiser le Lagrangien sans contrainte
#   Etape duale   : gradient projete pour les inegalites,
#                   gradient simple pour les egalites.
# ============================================================


# Documentation roxygen2 : ces lignes seront transformees
# automatiquement en fichier .Rd par devtools::document().
#' Algorithme d'Uzawa general (lineaire et non lineaire)
#'
#' Resout un probleme d'optimisation avec contraintes
#' d'inegalite et/ou d'egalite, lineaires ou non lineaires.
#'
#' @param f Fonction objectif f(x) renvoyant un scalaire.
#' @param x0 Point de depart (vecteur numeric).
#' @param g Fonction des contraintes d'inegalite g(x) <= 0
#'   renvoyant un vecteur (ou NULL si aucune).
#' @param h Fonction des contraintes d'egalite h(x) = 0
#'   renvoyant un vecteur (ou NULL si aucune).
#' @param rho Pas de mise a jour du multiplicateur (petit positif).
#' @param tol_var Tolerance sur la variation des multiplicateurs.
#' @param tol_feas Tolerance sur la violation des contraintes.
#' @param max_iter Nombre maximum d'iterations.
#' @param solver Methode passee a optim() : "BFGS", "L-BFGS-B", etc.
#' @param lower,upper Bornes optionnelles pour L-BFGS-B.
#' @param verbose Afficher la progression si TRUE.
#'
#' @return Une liste avec x, lambda, mu, iterations, converged.
#' @export
uzawa_general <- function(f, x0, g = NULL, h = NULL,
                          rho = 0.1,
                          tol_var = 1e-6, tol_feas = 1e-6,
                          max_iter = 1000,
                          solver = "BFGS",
                          lower = -Inf, upper = Inf,
                          verbose = FALSE) {

  # ----------------------------------------------------------
  # ETAPE 0 : determination du nombre de contraintes
  # ----------------------------------------------------------

  # Si g est NULL, m = 0 ; sinon on evalue g(x0) pour connaitre
  # le nombre de contraintes d'inegalite.
  m <- if (is.null(g)) 0L else length(g(x0))

  # Idem pour les contraintes d'egalite h.
  p <- if (is.null(h)) 0L else length(h(x0))

  # Si aucune contrainte n'est fournie, Uzawa n'a pas de sens :
  # on renvoie une erreur explicite.
  if (m == 0 && p == 0) {
    stop("Aucune contrainte fournie : utiliser optim() directement.")
  }

  # ----------------------------------------------------------
  # ETAPE 1 : initialisation
  # ----------------------------------------------------------

  # Multiplicateurs d'inegalite, initialises a zero.
  # Ils doivent rester >= 0 (contrainte duale).
  lambda <- rep(0, m)

  # Multiplicateurs d'egalite, initialises a zero.
  # Ils peuvent etre negatifs ou positifs (pas de projection).
  mu <- rep(0, p)

  # Point courant : on demarre au point fourni par l'utilisateur.
  x <- x0

  # ----------------------------------------------------------
  # ETAPE 2 : boucle principale de l'algorithme d'Uzawa
  # ----------------------------------------------------------
  for (k in seq_len(max_iter)) {

    # ---- Etape primale : minimisation du Lagrangien ----

    # On construit le Lagrangien comme une fonction de y
    # (variable muette), en figeant lambda et mu courants.
    lagrangien <- function(y) {
      # Terme objectif.
      val <- f(y)
      # Ajout des termes d'inegalite lambda_i * g_i(y).
      if (m > 0) val <- val + sum(lambda * g(y))
      # Ajout des termes d'egalite mu_j * h_j(y).
      if (p > 0) val <- val + sum(mu * h(y))
      # Valeur scalaire du Lagrangien au point y.
      val
    }

    # Resolution du probleme sans contrainte via optim().
    # On encapsule dans tryCatch pour eviter un crash si
    # l'optimiseur echoue.
    res <- tryCatch({
      # Cas avec bornes : on utilise L-BFGS-B.
      if (solver == "L-BFGS-B") {
        optim(par = x, fn = lagrangien, method = "L-BFGS-B",
              lower = lower, upper = upper,
              control = list(maxit = 500))
      } else {
        # Cas sans bornes : BFGS ou autre methode lisse.
        optim(par = x, fn = lagrangien, method = solver,
              control = list(maxit = 500))
      }
    }, error = function(e) NULL)

    # Si optim a echoue (res NULL), on s'arrete proprement.
    if (is.null(res)) {
      warning("Echec de l'etape primale a l'iteration ", k)
      return(list(x = x, lambda = lambda, mu = mu,
                  iterations = k, converged = FALSE))
    }

    # Nouveau point primal propose par optim().
    x_new <- res$par

    # ---- Etape duale : mise a jour des multiplicateurs ----

    # Par defaut, on conserve les anciens multiplicateurs.
    lambda_new <- lambda
    # S'il y a des inegalites : gradient projete sur R_+.
    # pmax(0, .) est la projection coordonnee par coordonnee.
    if (m > 0) {
      lambda_new <- pmax(0, lambda + rho * g(x_new))
    }

    # Pour les egalites : pas de projection (multiplicateur libre).
    mu_new <- mu
    if (p > 0) {
      mu_new <- mu + rho * h(x_new)
    }

    # ---- Criteres d'arret ----

    # Variation maximale des multiplicateurs entre deux iterations.
    delta_var <- max(
      if (m > 0) max(abs(lambda_new - lambda)) else 0,
      if (p > 0) max(abs(mu_new - mu))         else 0
    )

    # Violation maximale des contraintes au nouveau point.
    # Pour les inegalites : on ne penalise que g > 0.
    # Pour les egalites   : on prend la valeur absolue.
    viol_feas <- max(
      if (m > 0) max(pmax(0, g(x_new))) else 0,
      if (p > 0) max(abs(h(x_new)))     else 0
    )

    # Affichage de la progression si demande.
    if (verbose) {
      cat(sprintf("iter %4d | var = %.3e | viol = %.3e\n",
                  k, delta_var, viol_feas))
    }

    # Mise a jour effective des variables d'etat.
    lambda <- lambda_new
    mu     <- mu_new
    x      <- x_new

    # Test de convergence : les deux criteres doivent etre satisfaits.
    if (delta_var < tol_var && viol_feas < tol_feas) {
      return(list(x = x, lambda = lambda, mu = mu,
                  iterations = k, converged = TRUE,
                  violation = viol_feas))
    }
  }

  # Si on sort de la boucle sans convergence :
  warning("Nombre maximum d'iterations atteint.")
  list(x = x, lambda = lambda, mu = mu,
       iterations = max_iter, converged = FALSE,
       violation = viol_feas)
}
