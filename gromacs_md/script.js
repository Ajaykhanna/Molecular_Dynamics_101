/**
 * Toggles the visibility of a content element associated with a header.
 * Also toggles a 'collapsed' class on the header for styling (e.g., arrow direction).
 * @param {string} elementId - The ID of the content element to toggle.
 */
function toggleVisibility(elementId) {
    const element = document.getElementById(elementId);
    // Find the header element that calls this function with the specific elementId
    const header = document.querySelector(`[onclick="toggleVisibility('${elementId}')"]`);

    if (element && header) {
        // Use Tailwind's 'hidden' class to toggle visibility
        element.classList.toggle('hidden');
        // Toggle the 'collapsed' class on the header for CSS styling
        header.classList.toggle('collapsed');
    } else {
        console.warn(`Element with ID "${elementId}" or its corresponding header not found.`);
    }
}

// Initialize sections on page load.
// Hides elements linked to headers that initially have the 'collapsed' class.
document.addEventListener('DOMContentLoaded', function() {
    // Find all headers that should start collapsed
    const collapsedHeaders = document.querySelectorAll('.toggle-header.collapsed');

    collapsedHeaders.forEach(header => {
        // Extract the target element's ID from the onclick attribute
        const onclickAttr = header.getAttribute('onclick');
        if (onclickAttr) {
            // Basic regex to extract the ID between the quotes
            const match = onclickAttr.match(/toggleVisibility\('([^']+)'\)/);
            if (match && match[1]) {
                 const elementId = match[1];
                 const element = document.getElementById(elementId);
                 if (element) {
                     // Ensure the element starts hidden if the header is collapsed
                     element.classList.add('hidden');
                 } else {
                     console.warn(`Initial setup: Element with ID "${elementId}" not found for collapsed header.`);
                 }
            }
        }
    });
});
