document.addEventListener("DOMContentLoaded", function() {
  var images = document.querySelectorAll("body img");
  images.forEach(function(img) {
    img.style.cursor = "pointer";
    img.onclick = function() {
      showImageInModal(img);
    };
  });
});

function showImageInModal(imgElement) {
  var imgSrc = imgElement.getAttribute("src");
  var modal = document.createElement("div");
  modal.style.position = "fixed";
  modal.style.top = "0";
  modal.style.left = "0";
  modal.style.width = "100%";
  modal.style.height = "100%";
  modal.style.backgroundColor = "rgba(0,0,0,0.8)";
  modal.style.display = "flex";
  modal.style.justifyContent = "center";
  modal.style.alignItems = "center";
  modal.style.zIndex = "1000";

  var fullImage = document.createElement("img");
  fullImage.src = imgSrc;
  fullImage.style.maxWidth = "100%";
  fullImage.style.maxHeight = "100%";
  modal.appendChild(fullImage);

  modal.onclick = function() {
    document.body.removeChild(modal);
  };

  document.body.appendChild(modal);
}
